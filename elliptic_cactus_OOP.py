import os
import shutil
import json
import math
import requests
from datetime import datetime, timedelta

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy

from PIL import Image
from skimage.measure import label, regionprops
from matplotlib.patches import Ellipse
from skimage import exposure
from skimage.filters import meijering
from skimage.morphology import disk, closing

# Assuming this is your local import
from extractDataSrc.Lasco import Lasco


class CMEDetector:
    def __init__(self, init_date: str, download_images: bool = False, base_folder: str = './data_processed/lasco/c3/'):
        """
        Initializes the CME Detector and runs the entire algorithm.
        
        :param init_date: The onset event date (e.g., '2023-01-01 12:00:00')
        :param download_images: Boolean flag to download new images or use existing ones.
        :param base_folder: The directory where images and metadata are stored.
        """
        self.init_date = init_date
        self.download_images = download_images
        self.base_folder = base_folder
        
        # Constants
        self.ARCSEC_TO_KM = 725.0    # Approximate conversion at 1 AU
        self.V_MIN_KM_S = 250.0      # Minimum plausible CME speed
        self.V_MAX_KM_S = 2000.0     # Maximum plausible CME speed
        self.CENTER_X = 512
        self.CENTER_Y = 512
        
        # Placeholders for data
        self.images = []
        self.metadata = {}
        self.arcsec_per_pixel = None
        self.rsun_pixels = None
        self.phi = None
        
        # Result containers
        self.detections = pd.DataFrame()
        self.cluster_dict = {}
        self.best_detection_dict = {}
        self.best_cluster_QS = 0.0

        # ---------------------------------------------------------
        # Execute the algorithm pipeline on instantiation
        # ---------------------------------------------------------
        self._run_pipeline()

    def _get_date_range(self, date_str):
        dt = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S')
        start = (dt - timedelta(hours=2)).replace(minute=0, second=0)
        end = dt + timedelta(hours=12)
        return start.strftime('%Y-%m-%d %H:%M:%S'), end.strftime('%Y-%m-%d %H:%M:%S')

    def _download_metadata(self, time_str):
        original_datetime = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
        date_str = original_datetime.strftime("%Y-%m-%dT%H:%M:%S.000Z")
        url = f"https://api.helioviewer.org/v2/getClosestImage/?date={date_str}&sourceId=5"

        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            os.makedirs(self.base_folder, exist_ok=True)
            with open(os.path.join(self.base_folder, "images_metadata.json"), 'w') as f:
                json.dump(data, f, indent=4)
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data from Helioviewer API: {e}")
        except IOError as e:
            print(f"Error saving metadata file: {e}")

    def _prepare_environment(self):
        if self.download_images:
            print('Downloading coronagraph images -----------------')
            if os.path.exists(self.base_folder):
                for filename in os.listdir(self.base_folder):
                    file_path = os.path.join(self.base_folder, filename)
                    try:
                        if os.path.isfile(file_path) or os.path.islink(file_path):
                            os.unlink(file_path)
                        elif os.path.isdir(file_path):
                            shutil.rmtree(file_path)
                    except Exception as e:
                        print(f'Failed to delete {file_path}. Reason: {e}')
            
            start_date, stop_date = self._get_date_range(self.init_date)
            self._download_metadata(start_date)
            lasco = Lasco(start_date, stop_date)
            lasco.extract_data("c3")

    def _load_data(self):
        for im in sorted(os.listdir(self.base_folder)):
            if im.endswith(".png"):
                img = Image.open(os.path.join(self.base_folder, im))
                self.images.append(np.array(img.convert('L')))
        
        json_file = os.path.join(self.base_folder, 'images_metadata.json')
        with open(json_file, 'r') as f:
            self.metadata = json.load(f)
            
        self.arcsec_per_pixel = self.metadata['scale']
        self.rsun_pixels = self.metadata['rsun']

    def _calculate_pylon(self):
        if not self.images:
            raise ValueError("No images found to calculate pylon.")
            
        first_image_data = self.images[0]
        height, width = first_image_data.shape

        y_indices, x_indices = np.indices((height, width))
        x_rel_pixels = x_indices - self.CENTER_X
        y_rel_pixels = self.CENTER_Y - y_indices

        x_arcsec = x_rel_pixels * self.arcsec_per_pixel
        y_arcsec = y_rel_pixels * self.arcsec_per_pixel
        r_arcsec = np.sqrt(x_arcsec**2 + y_arcsec**2)
        theta_arcsec = np.mod(np.arctan2(y_arcsec, x_arcsec), 2 * np.pi)

        best_pylon_cover = None
        for pylon_angle in np.arange(0.01, 2, 0.01):
            pylon_mask = (theta_arcsec < pylon_angle * np.pi) & \
                         (theta_arcsec > (pylon_angle - 7/100) * np.pi) & \
                         (r_arcsec > 4.7 * self.rsun_pixels * self.arcsec_per_pixel) & \
                         (r_arcsec < 512 * self.arcsec_per_pixel)

            pylon = first_image_data[pylon_mask]
            if len(pylon) == 0:
                continue
                
            match_count = np.sum(pylon < 5)
            total_count = len(pylon)

            if best_pylon_cover is None or match_count / total_count > best_pylon_cover:
                best_pylon_cover = match_count / total_count
                self.phi = pylon_angle

    def _prepare_image(self, image):
        height, width = image.shape
        y_indices, x_indices = np.indices((height, width))

        x_rel_pixels = x_indices - self.CENTER_X
        y_rel_pixels = y_indices - self.CENTER_Y

        x_arcsec = x_rel_pixels * self.arcsec_per_pixel
        y_arcsec = y_rel_pixels * self.arcsec_per_pixel

        r_arcsec = np.sqrt(x_arcsec**2 + y_arcsec**2)
        theta_arcsec = np.arctan2(y_arcsec, x_arcsec)
        
        # Mask out sun, occulter, edges
        mask = (r_arcsec < 4.7 * self.rsun_pixels * self.arcsec_per_pixel) | \
               (r_arcsec > 512 * self.arcsec_per_pixel) | \
               ((theta_arcsec < self.phi * np.pi) & (theta_arcsec > (self.phi - 7/100) * np.pi))

        image_data_masked = image.copy().astype(np.float32)
        image_data_masked[mask] = np.nan

        return image_data_masked, r_arcsec, theta_arcsec

    def _unroll_to_polar(self, image_data, r_arcsec_array, theta_arcsec_array, num_radii=512, num_angles=360):
        R_max = r_arcsec_array.max()
        R_new = np.linspace(r_arcsec_array.min(), R_max, num_radii)
        Theta_new = np.linspace(0, 2 * np.pi, num_angles, endpoint=False)
        R_grid, Theta_grid = np.meshgrid(R_new, Theta_new)
        
        indices_x = self.CENTER_X + R_grid * np.cos(Theta_grid) / self.arcsec_per_pixel
        indices_y = self.CENTER_Y + R_grid * np.sin(Theta_grid) / self.arcsec_per_pixel

        coords = np.vstack([indices_y.ravel(), indices_x.ravel()])

        polar_map = scipy.ndimage.map_coordinates(
            input=image_data, 
            coordinates=coords, 
            order=1,
            cval=np.nan 
        ).reshape(num_angles, num_radii)
        
        return polar_map.T

    def _create_jmap(self, diff_ims, angle):
        j_map_slices = [polar_map[:, angle] for polar_map in diff_ims[::-1]]
        return np.stack(j_map_slices, axis=0)

    def _cme_blob_detection(self, j_map, percentil, disk_size, verbose=False):
        j_map_clean = np.nan_to_num(j_map, nan=0.0, posinf=0.0, neginf=0.0)
        p1, p99 = np.percentile(j_map_clean, (1, 99))
        j_map_clipped = np.clip(j_map_clean, p1, p99)
        j_map_norm = (j_map_clipped - j_map_clipped.min()) / (j_map_clipped.max() - j_map_clipped.min())
        j_map_eq = exposure.equalize_adapthist(j_map_norm, clip_limit=0.03)

        j_map_blobs = meijering(j_map_eq, sigmas=range(1, 4), black_ridges=False)
        threshold = np.percentile(j_map_blobs, percentil) 
        binary_j_map = closing(j_map_blobs > threshold, disk(disk_size))
        label_image = label(binary_j_map)

        if verbose:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.imshow(binary_j_map, cmap='gray', origin='upper')

        cme_detection = []

        for props in regionprops(label_image):
            if props.area < 40:
                continue

            y0, x0 = props.centroid
            major_axis = props.major_axis_length
            minor_axis = props.minor_axis_length
            orientation_rad = props.orientation
            angle_rad = orientation_rad + (np.pi / 2)
            
            x_coords = props.coords[:, 1]
            y_coords = props.coords[:, 0]
            unique_radii = np.unique(x_coords)

            min_y_edge = [np.min(y_coords[x_coords == r]) for r in unique_radii]
            max_y_edge = [np.max(y_coords[x_coords == r]) for r in unique_radii]

            slope_min = np.polyfit(unique_radii, min_y_edge, 1)[0]
            slope_max = np.polyfit(unique_radii, max_y_edge, 1)[0]
            slope = min(abs(slope_min), abs(slope_max))
            
            track_length_x = np.max(x_coords) - np.min(x_coords)
            if track_length_x < 30 or slope <= 0.02:
                continue

            R_arcsec_max = (1024 / 2) * self.arcsec_per_pixel
            R_arcsec_min = 3.7 * self.rsun_pixels * self.arcsec_per_pixel
            dr_arcsec_per_pixel = (R_arcsec_max - R_arcsec_min) / 512 
            dt_seconds = 12 * 60.0 

            velocity_arcsec_per_sec = abs(1.0 / slope) * (dr_arcsec_per_pixel / dt_seconds)
            velocity_km_s = velocity_arcsec_per_sec * self.ARCSEC_TO_KM
            t_onset_idx = abs(y0 + slope * x0)
            
            if self.V_MIN_KM_S < velocity_km_s < self.V_MAX_KM_S:
                cme_detection.append({
                    'x0': x0, 'y0': y0,
                    'angle_deg': math.degrees(angle_rad),
                    'major_axis': major_axis, 'minor_axis': minor_axis,
                    'slope': slope, 'velocity_km_s': velocity_km_s,
                    't_onset_idx': math.ceil(t_onset_idx)
                })

        return cme_detection

    def _all_angles_detections_blobs(self, images):
        rows = []
        for angle in range(0, 360, 5):
            j_map = self._create_jmap(images, angle)
            cme = self._cme_blob_detection(j_map, percentil=97, disk_size=2, verbose=False)
            for r in cme:
                r['angle'] = angle
                rows.append(r)
        return pd.DataFrame(rows)

    def _fill_factor(self, cluster):
        c1_angles = cluster['angle'].sort_values().values
        if len(c1_angles) <= 1:
            return 0.0, 0.0
            
        gaps = np.diff(c1_angles)
        wrap_gap = (360 - c1_angles[-1]) + c1_angles[0]
        all_gaps = np.append(gaps, wrap_gap)
        c1_AW = 360 - np.max(all_gaps)
        
        if c1_AW == 0:
            return 0.0, 0.0
            
        FF_c1 = (len(c1_angles) / ((c1_AW / 5) + 1)) * 100
        return min(100.0, FF_c1), c1_AW

    def _quality_score(self, flat_clusters):
        cluster_dict = {}
        # Fetching images specifically sorted
        img_files = sorted([im for im in os.listdir(self.base_folder) if im.endswith('.png')])

        if flat_clusters.empty:
            return {}

        for cluster_id in flat_clusters['cme_cluster_id'].unique():
            cluster = flat_clusters[flat_clusters['cme_cluster_id'] == cluster_id]
            cluster = cluster.sort_values('angle')
            smoothed_velocities = cluster['velocity_km_s'].rolling(window=3, center=True, min_periods=1).median()
            true_cme_velocity = smoothed_velocities.max()
            
            v = cluster['velocity_km_s']
            t = cluster['t_onset_idx']

            cv_t = t.std() / t.mean() if t.mean() != 0 else 1.0
            cv_v = v.std() / v.mean() if v.mean() != 0 else 1.0

            FF, AW = self._fill_factor(cluster)
            norm_AW = AW / 360.0
            norm_FF = FF / 100.0
            norm_cv_v = max(0.0, 1.0 - cv_v)
            norm_cv_t = max(0.0, 1.0 - cv_t)

            qs = (0.4 * norm_AW) + (0.3 * norm_FF) + (0.15 * norm_cv_v) + (0.15 * norm_cv_t)

            ot = math.floor(t.mean().item())
            iot = len(self.images) - ot
            
            try:
                onset_image_date, onset_image_datetime, _, _ = img_files[iot].split('_')
            except (IndexError, ValueError):
                onset_image_date, onset_image_datetime = "Unknown", "Unknown"

            cluster_dict[cluster_id.item()] = {
                "zscore_velocity": cv_v.item() if not pd.isna(cv_v.item()) else 0.0,
                "zscore_onset_time": cv_t.item() if not pd.isna(cv_t.item()) else 0.0,
                "velocity_km_s": true_cme_velocity.item(),
                "onset_time_idx": ot,
                'onset_time_inverse_idx': iot,
                'onset_date': onset_image_date,
                'onset_datetime': onset_image_datetime,
                "angular_width": AW.item(),
                'fill_factor': FF.item(),
                'THETA': qs.item()
            }

        return cluster_dict

    def _run_pipeline(self):
        """
        Executes the main pipeline. 
        All operations originally at the bottom of the script are here.
        """
        self._prepare_environment()
        self._load_data()
        self._calculate_pylon()

        print("Preparing images and creating polar maps...")
        prepared_images, prepared_r_arcsecs, prepared_theta_arcsecs = [], [], []
        
        for image in self.images:
            im, r, t = self._prepare_image(image)
            prepared_images.append(im)
            prepared_r_arcsecs.append(r)
            prepared_theta_arcsecs.append(t)

        polar_maps = [
            self._unroll_to_polar(img, r, t) 
            for img, r, t in zip(prepared_images, prepared_r_arcsecs, prepared_theta_arcsecs)
        ]

        print("Calculating running differences...")
        running_difference = [
            np.clip(polar_maps[i] - polar_maps[i-1], 0, 255) 
            for i in range(1, len(polar_maps))
        ]

        print("Running blob detections...")
        detections = self._all_angles_detections_blobs(running_difference)
        
        if detections.empty:
            print("No CMEs detected.")
            return

        detections = detections.sort_values(by='t_onset_idx').reset_index(drop=True)
        detections['onset_diff'] = detections['t_onset_idx'].diff()
        detections.fillna(0, inplace=True)
        detections['cme_cluster_id'] = (detections['onset_diff'] > 1).cumsum()

        counts = detections['cme_cluster_id'].value_counts()
        valid_clusters = counts[counts > 10].index
        self.detections = detections[detections['cme_cluster_id'].isin(valid_clusters)]

        flat_clusters = pd.DataFrame()
        for c in self.detections['cme_cluster_id'].unique():
            cc = self.detections[self.detections['cme_cluster_id'] == c].groupby('angle').agg({
                't_onset_idx': 'median',
                'velocity_km_s': 'median',
                'cme_cluster_id': 'min'
            }).reset_index()
            flat_clusters = pd.concat([flat_clusters, cc])

        print("Calculating quality scores...")
        self.cluster_dict = self._quality_score(flat_clusters)
        
        print("\n--- Process Complete ---")
        print(json.dumps(self.cluster_dict, indent=4))

        if not self.cluster_dict:
            for d in self.cluster_dict:
                if self.cluster_dict[d]['THETA'] > self.best_cluster_QS:
                    self.best_cluster_QS = self.cluster_dict[d]['THETA']
                    self.best_detection_dict = self.cluster_dict[d]

# ==========================================
# How to use the class
# ==========================================
if __name__ == "__main__":
    # Create the object, pass the date, and tell it whether to download images
    # The pipeline executes automatically upon creation.

    want_download = input("Do you want to download new images? (y/n): ").strip().lower() == 'y'
    date = None
    if want_download:
        date = input("Enter the onset event date (e.g., '2024-05-09 09:24:00'): ").strip()

    detector = CMEDetector(init_date=date, download_images=want_download)

    print("Best Dict: ")
    print(json.dumps(
        detector.best_detection_dict, indent=4
    ))
    
    # You can access your final data directly from the object instance
    # print(analyzer.detections.head())
    # print(analyzer.cluster_dict)