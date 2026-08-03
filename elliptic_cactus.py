import os, shutil

import requests
from extractDataSrc.Lasco import Lasco
from datetime import datetime, timedelta

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import scipy

from PIL import Image
import json


from skimage.measure import label, regionprops
import math
from matplotlib.patches import Ellipse

from skimage import exposure
from skimage.filters import meijering, frangi, sato
from skimage.morphology import disk, closing, footprint_rectangle

from sklearn.cluster import DBSCAN
from sklearn.linear_model import RANSACRegressor
from scipy.interpolate import interp1d

import warnings




# ─────────────────────────────────────────────
#  Constants (physical, not path-dependent)
# ─────────────────────────────────────────────
ARCSEC_TO_KM = 725.0    # Approximate conversion at 1 AU
V_MIN_KM_S   = 250.0    # Minimum plausible CME speed
V_MAX_KM_S   = 2000.0   # Maximum plausible CME speed


# ─────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────

def get_date_range(date_str):
    dt = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S')
    start = (dt - timedelta(hours=2)).replace(minute=0, second=0)
    end   = dt + timedelta(hours=12)
    return (
        start.strftime('%Y-%m-%d %H:%M:%S'),
        end.strftime('%Y-%m-%d %H:%M:%S')
    )


def load_images(dir_path):
    """Load all PNG images from dir_path.
    Expects a JSON metadata file called images_metadata.json in the same directory.
    """
    images = []
    for im in sorted(os.listdir(dir_path)):
        if im.endswith(".png"):
            img = Image.open(os.path.join(dir_path, im))
            images.append(np.array(img.convert('L')))

    # Use the first PNG as the reference image for shape
    first_png = next(f for f in sorted(os.listdir(dir_path)) if f.endswith('.png'))
    img = Image.open(os.path.join(dir_path, first_png))
    image_data = np.array(img.convert('L'))

    # Metadata lives alongside the images
    json_file = os.path.join(dir_path, 'images_metadata.json')
    with open(json_file, 'r') as f:
        metadata = json.load(f)

    return images, image_data, metadata


def download_metadata(time_str, save_dir):
    original_datetime = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
    date_str = original_datetime.strftime("%Y-%m-%dT%H:%M:%S.000Z")
    url = f"https://api.helioviewer.org/v2/getClosestImage/?date={date_str}&sourceId=5"

    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        out_path = os.path.join(save_dir, "images_metadata.json")
        with open(out_path, 'w') as f:
            json.dump(data, f, indent=4)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from Helioviewer API: {e}")
    except IOError as e:
        print(f"Error saving metadata: {e}")


# ─────────────────────────────────────────────
#  Core pipeline – everything that used to be
#  at module level is now inside this function.
# ─────────────────────────────────────────────

def run_pipeline(dir_path):
    """Run the full CME detection pipeline for images in dir_path.

    Args:
        dir_path: Path to a directory containing PNGs + images_metadata.json.

    Returns:
        dict: cluster_dict produced by quality_score().
    """
    print(f"\n{'='*60}")
    print(f"Running pipeline for: {dir_path}")
    print(f"{'='*60}")

    images, image_data, metadata = load_images(dir_path)

    # ── Coordinate grids (depend on metadata, so computed per-run) ──
    center_x = 512
    center_y = 512

    arcsec_per_pixel = metadata['scale']
    rsun_pixels      = metadata['rsun']

    height, width = image_data.shape
    y_indices, x_indices = np.indices((height, width))

    x_rel_pixels = x_indices - center_x
    y_rel_pixels = center_y - y_indices

    x_arcsec = x_rel_pixels * arcsec_per_pixel
    y_arcsec = y_rel_pixels * arcsec_per_pixel

    r_arcsec     = np.sqrt(x_arcsec**2 + y_arcsec**2)
    theta_arcsec = np.arctan2(y_arcsec, x_arcsec)
    theta_arcsec = np.mod(theta_arcsec, 2 * np.pi)

    R_MIN_ARCSEC = 4.7 * rsun_pixels * arcsec_per_pixel
    R_MAX_ARCSEC = 512 * arcsec_per_pixel

    # ── Pylon angle detection ──
    best_pylon_cover = None
    phi = None
    for pylon_angle in np.arange(0.01, 2, 0.01):
        pylon_mask = (
            (theta_arcsec < pylon_angle * np.pi) &
            (theta_arcsec > (pylon_angle - 7/100) * np.pi) &
            (r_arcsec > 4.7 * rsun_pixels * arcsec_per_pixel) &
            (r_arcsec < 512 * arcsec_per_pixel)
        )
        pylon       = image_data[pylon_mask]
        match_count = np.sum(pylon < 5)
        total_count = len(pylon)

        if best_pylon_cover is None or match_count / total_count > best_pylon_cover:
            best_pylon_cover = match_count / total_count
            phi = pylon_angle

    # ── Inner functions that close over the per-run variables ──

    def image_masking(img):
        mask = (
            (r_arcsec < 4.7 * rsun_pixels * arcsec_per_pixel) |
            (r_arcsec > 512 * arcsec_per_pixel) |
            ((theta_arcsec < phi * np.pi) & (theta_arcsec > (phi - 7/100) * np.pi))
        )
        img_masked = img.copy().astype(np.float32)
        img_masked[mask] = np.nan
        return img_masked

    def prepare_image(image):
        h, w = image_data.shape
        yi, xi = np.indices((h, w))

        xr = (xi - center_x) * arcsec_per_pixel
        yr = (yi - center_y) * arcsec_per_pixel

        r_local     = np.sqrt(xr**2 + yr**2)
        theta_local = np.arctan2(yr, xr)

        mask = (
            (r_local < 4.7 * rsun_pixels * arcsec_per_pixel) |
            (r_local > 512 * arcsec_per_pixel) |
            ((theta_local < phi * np.pi) & (theta_local > (phi - 7/100) * np.pi))
        )
        img_masked = image.copy().astype(np.float32)
        img_masked[mask] = np.nan
        return img_masked, r_local, theta_local

    def unroll_to_polar(img, r_arcsec_array, theta_arcsec_array,
                        num_radii=512, num_angles=360):
        R_max   = r_arcsec_array.max()
        R_new = np.linspace(R_MIN_ARCSEC, R_MAX_ARCSEC, num_radii)
        Theta_new = np.linspace(0, 2 * np.pi, num_angles, endpoint=False)
        R_grid, Theta_grid = np.meshgrid(R_new, Theta_new)

        occulter_cx = 512.0
        occulter_cy = 512.0
        indices_x = occulter_cx + R_grid * np.cos(Theta_grid) / arcsec_per_pixel
        indices_y = occulter_cy + R_grid * np.sin(Theta_grid) / arcsec_per_pixel

        coords = np.vstack([indices_y.ravel(), indices_x.ravel()])
        polar_map = scipy.ndimage.map_coordinates(
            input=img,
            coordinates=coords,
            order=1,
            cval=np.nan
        ).reshape(num_angles, num_radii)

        return polar_map.T

    def create_jmap(diff_ims, center_angle, angular_width=5):
        j_map_slices = []
        half_width = angular_width // 2
        max_angle  = diff_ims[0].shape[1]
        start_angle_idx = max(0, center_angle - half_width)
        end_angle_idx   = min(max_angle, center_angle + half_width + 1)

        for polar_map in diff_ims[::-1]:
            swath = polar_map[:, start_angle_idx:end_angle_idx]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                clean_radial_slice = np.nanmean(swath, axis=1)
            j_map_slices.append(clean_radial_slice)

        return np.stack(j_map_slices, axis=0)

    def resample_jmap_to_constant_cadence(j_map, cum_times_seconds, target_cadence_mins=12):
        """
        Interpolates a J-map with variable time gaps to a constant time step.
        
        j_map: The original 2D array [Time, Space]
        cum_times_seconds: 1D array of exact seconds for each row in j_map
        target_cadence_mins: The constant time step you want (e.g., 12 mins)
        """
        # The original time axis in minutes
        original_times_mins = (cum_times_seconds / 60.0)
        
        # Create a new, perfectly even time axis (every 12 minutes)
        max_time = original_times_mins[-1]
        new_times_mins = np.arange(0, max_time + target_cadence_mins, target_cadence_mins)
        
        # Interpolate the J-map along the time axis (axis 0)
        # interp1d creates a mathematical function from your data
        print("O"*50 + f"{len(original_times_mins)}, {j_map.shape}")
        interpolator = interp1d(original_times_mins, j_map, axis=0, kind='linear', fill_value="extrapolate")
        
        # Generate the new J-map using the perfectly even time steps
        resampled_jmap = interpolator(new_times_mins)
        
        return resampled_jmap

    def cme_blob_detection(j_map, percentil, disk_size, verbose=False):
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
            major_axis = props.axis_major_length
            minor_axis = props.axis_minor_length
            orientation_rad = props.orientation
            angle_rad = orientation_rad + (np.pi / 2)
            
            x_coords = props.coords[:, 1]
            y_coords = props.coords[:, 0]
            unique_radii = np.unique(x_coords)

            min_y_edge = [np.min(y_coords[x_coords == r]) for r in unique_radii]
            max_y_edge = [np.max(y_coords[x_coords == r]) for r in unique_radii]


            # slope_min = np.polyfit(unique_radii, min_y_edge, 1)[0]
            # slope_max = np.polyfit(unique_radii, max_y_edge, 1)[0]
            # slope = min(abs(slope_min), abs(slope_max))

            # -------------------------------------------------------------
            # THE RANSAC FIX (Restored)
            # -------------------------------------------------------------
            # Reshape X for scikit-learn
            X = unique_radii.reshape(-1, 1)
            
            try:
                # residual_threshold=1.0 forces it to ignore the "fat" from disk()
                ransac_min = RANSACRegressor(residual_threshold=1.0, random_state=42)
                ransac_min.fit(X, min_y_edge)
                slope_min = ransac_min.estimator_.coef_[0]

                ransac_max = RANSACRegressor(residual_threshold=1.0, random_state=42)
                ransac_max.fit(X, max_y_edge)
                slope_max = ransac_max.estimator_.coef_[0]
                
                slope = min(abs(slope_min), abs(slope_max))
                
            except ValueError:
                # Fallback to polyfit if the blob is too small for RANSAC
                slope_min = np.polyfit(unique_radii, min_y_edge, 1)[0]
                slope_max = np.polyfit(unique_radii, max_y_edge, 1)[0]
                slope = min(abs(slope_min), abs(slope_max))
            # ------------------------------------------------------------
            
            track_length_x = np.max(x_coords) - np.min(x_coords)
            if track_length_x < 30 or slope <= 0.02:
                continue

            # R_arcsec_max = (1024 / 2) * arcsec_per_pixel
            # R_arcsec_min = 3.7 * rsun_pixels * arcsec_per_pixel
            # dr_arcsec_per_pixel = (R_arcsec_max - R_arcsec_min) / 512 
            dr_arcsec_per_pixel = (R_MAX_ARCSEC - R_MIN_ARCSEC) / 512
            dt_seconds = 12 * 60.0 

            velocity_arcsec_per_sec = abs(1.0 / slope) * (dr_arcsec_per_pixel / dt_seconds)
            velocity_km_s = velocity_arcsec_per_sec * ARCSEC_TO_KM
            t_onset_idx = abs(y0 + slope * x0)
            
            if V_MIN_KM_S < velocity_km_s < V_MAX_KM_S:
                cme_detection.append({
                    'x0': x0, 'y0': y0,
                    'angle_deg': math.degrees(angle_rad),
                    'major_axis': major_axis, 'minor_axis': minor_axis,
                    'slope': slope, 'velocity_km_s': velocity_km_s,
                    't_onset_idx': math.ceil(t_onset_idx)
                })

        return cme_detection

    def all_angles_detections_blobs(images_in, cum_times):
        rows = []
        for angle in range(0, 360, 5):
            j_map = create_jmap(images_in, angle)
            
            total_duration_sec = cum_times[-1] - cum_times[0]
            expected_intervals = total_duration_sec / 720
            actual_intervals = len(cum_times) - 1
            print(50*'-' + f"{actual_intervals}, {j_map.shape}")
            cadence_score = actual_intervals / expected_intervals
            
            if cadence_score < 0.7:
                # FIX 1: Drop the first timestamp to match the running difference length (N-1)
                cum_times_diff = cum_times[1:]
                
                original_jmap_len = j_map.shape[0]
                
                # Pass the sliced time array
                j_map_resampled = resample_jmap_to_constant_cadence(j_map=j_map, cum_times_seconds=cum_times_diff)
                resampled_jmap_len = j_map_resampled.shape[0]
                
                cme = cme_blob_detection(j_map_resampled, percentil=97, disk_size=2, verbose=False)
                
                # FIX 2: Scale the onset index back to original J-map space
                scale_factor = original_jmap_len / resampled_jmap_len
                
                for r in cme:
                    r['angle'] = angle
                    # Convert the resampled index back to normal index so the dates work
                    r['t_onset_idx'] = r['t_onset_idx'] * scale_factor
                    rows.append(r)
            else:
                # Cadence is good, run normally
                cme = cme_blob_detection(j_map, percentil=97, disk_size=2, verbose=False)
                for r in cme:
                    r['angle'] = angle
                    rows.append(r)
                    
        return pd.DataFrame(rows)

    def fill_factor(cluster):
        c1_angles = cluster['angle'].sort_values().values
        if len(c1_angles) <= 1:
            return 0.0, 0.0
        gaps     = np.diff(c1_angles)
        wrap_gap = (360 - c1_angles[-1]) + c1_angles[0]
        all_gaps = np.append(gaps, wrap_gap)
        c1_AW    = 360 - np.max(all_gaps)
        if c1_AW == 0:
            return 0.0, 0.0
        FF_c1 = (len(c1_angles) / ((c1_AW / 5) + 1)) * 100
        return min(100.0, FF_c1), c1_AW
    
    def estimate_angular_width(diff_ims, percentile_threshold=95):
        """
        Estimates the angular width of the CME from polar difference images.
        NaN-safe to handle occulter and pylon masks.
        """
        cube = np.stack(diff_ims, axis=0)
        
        # 1. Use nanmax to ignore the NaNs from the pylon and edges!
        # We use warnings.catch_warnings to suppress the warning Numpy gives 
        # when it encounters a slice that is *entirely* NaNs.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            max_intensity_profile = np.nanmax(cube, axis=(0, 1))
            
        # 2. Use nanpercentile for the threshold
        threshold = np.nanpercentile(max_intensity_profile, percentile_threshold)
        
        active_angles = np.where(max_intensity_profile > threshold)[0]
        
        if len(active_angles) <= 1:
            return 0 # No significant CME found
            
        num_angles = diff_ims[0].shape[1] # Usually 360
        
        # 3. Calculate internal gaps
        diffs = np.diff(active_angles)
        
        # 4. Calculate the gap that wraps across the 360 -> 0 line
        wrap_gap = (num_angles - active_angles[-1]) + active_angles[0]
        
        # 5. Combine all gaps to find the true largest empty space
        all_gaps = np.append(diffs, wrap_gap)
        max_gap = np.max(all_gaps)
        
        # If the largest empty gap is greater than 90 degrees, it's a partial CME.
        if max_gap > (num_angles / 4): 
            width_pixels = num_angles - max_gap
        else:
            width_pixels = num_angles
            
        degrees_per_pixel = 360.0 / num_angles
        angular_width_deg = width_pixels * degrees_per_pixel
        
        return angular_width_deg
    
    def correct_projection_effect(velocity_pos_km_s, angular_width_deg):
        """
        Reduces the apparent Plane-of-Sky speed for wide Halo CMEs 
        to approximate the true 3D radial speed.
        """
        if angular_width_deg <= 120:
            return velocity_pos_km_s

        halo_factor = min((angular_width_deg - 120) / 240.0, 1.0)
        max_penalty = 0.15  # Up to 15% reduction for a full 360 Halo
        correction_multiplier = 1.0 - (halo_factor * max_penalty)

        return velocity_pos_km_s * correction_multiplier

    def quality_score(flat_clusters, diff_ims):
        cluster_dict = {}

        if flat_clusters.empty:
            return { 
                '0': {
                    "zscore_velocity":    0,
                    "zscore_onset_time":  0,
                    "velocity_km_s":      0,
                    "onset_time_idx":     0,
                    'onset_time_inverse_idx': 0,
                    'onset_date':         0,
                    'onset_datetime':     0,
                    "angular_width":      0,
                    'fill_factor':        0,
                    'THETA':             -1
                }
            }

        # Sorted list of PNG filenames for onset-time look-up
        png_files = sorted(f for f in os.listdir(dir_path) if f.endswith('.png'))

        for cluster_id in flat_clusters['cme_cluster_id'].unique():
            cluster = flat_clusters[flat_clusters['cme_cluster_id'] == cluster_id]

            cluster = cluster.sort_values('angle')  
            smoothed_velocities = cluster['velocity_km_s'].rolling(window=3, center=True, min_periods=1).median()
            true_cme_velocity = smoothed_velocities.max()

            FF, __AW = fill_factor(cluster)
            AW = estimate_angular_width(diff_ims)
            true_cme_velocity_corrected = correct_projection_effect(true_cme_velocity, __AW.item())

            v = cluster['velocity_km_s']
            t = cluster['t_onset_idx']

            cv_t = t.std() / t.mean() if t.mean() != 0 else 1.0
            cv_v = v.std() / v.mean() if v.mean() != 0 else 1.0


            norm_AW   = __AW / 360.0
            norm_FF   = FF / 100.0
            norm_cv_v = max(0.0, 1.0 - cv_v)
            norm_cv_t = max(0.0, 1.0 - cv_t)

            qs = (0.4 * norm_AW) + (0.3 * norm_FF) + (0.15 * norm_cv_v) + (0.15 * norm_cv_t)

            ot  = math.floor(t.mean().item())
            iot = len(images) - ot

            onset_image_date, onset_image_datetime, _, _ = png_files[iot].split('_')

            cluster_dict[cluster_id.item()] = {
                "zscore_velocity":        cv_v.item(),
                "zscore_onset_time":      cv_t.item(),
                "velocity_km_s":          true_cme_velocity_corrected.item(),
                "onset_time_idx":         ot,
                'onset_time_inverse_idx': iot,
                'onset_date':             onset_image_date,
                'onset_datetime':         onset_image_datetime,
                "angular_width":          __AW.item(), ## in the case __AW use __AW.item()
                'fill_factor':            FF,
                'THETA':                  qs.item()
            }

        return cluster_dict

    # ── Image preparation ──
    print("-"*50)
    print("PREPARING IMAGES")
    print("-"*50)
    prepared_images      = []
    prepared_r_arcsecs   = []
    prepared_theta_arcsecs = []

    for image in images:
        im, r, t = prepare_image(image)
        prepared_images.append(im)
        prepared_r_arcsecs.append(r)
        prepared_theta_arcsecs.append(t)

    # ── Cumulative times ──
    files = sorted(f for f in os.listdir(dir_path) if f.endswith('.png'))
    time_strings = []
    for im in files:
        parts         = im.split('_')
        full_time_str = f"{parts[0]}_{parts[1]}"
        time_strings.append(full_time_str)

    times    = pd.to_datetime(time_strings, format='%Y%m%d_%H%M')
    t0       = times[0]
    cum_times = np.array([(t - t0).total_seconds() for t in times])

    # ── Polar transform + running difference ──
    polar_maps = [
        unroll_to_polar(img, r, t)
        for img, r, t in zip(prepared_images, prepared_r_arcsecs, prepared_theta_arcsecs)
    ]

    running_difference = [
        np.clip(polar_maps[i] - polar_maps[i-1], 0, 255)
        for i in range(1, len(polar_maps))
    ]

    # ── Detection ──
    print("-"*50)
    print("DETECTION")
    print("-"*50)
    detections = all_angles_detections_blobs(running_difference, cum_times=cum_times)
    detections = detections.sort_values(by='t_onset_idx').reset_index(drop=True)

    detections['onset_diff']    = detections['t_onset_idx'].diff()
    detections.fillna(0, inplace=True)
    detections['cme_cluster_id'] = (detections['onset_diff'] > 1).cumsum()

    counts         = detections['cme_cluster_id'].value_counts()
    valid_clusters = counts[counts > 10].index
    detections     = detections[detections['cme_cluster_id'].isin(valid_clusters)]

    flat_clusters = pd.DataFrame()
    for c in detections['cme_cluster_id'].unique():
        cc = detections[detections['cme_cluster_id'] == c].groupby('angle').agg({
            't_onset_idx':   'median',
            'velocity_km_s': 'median',
            'cme_cluster_id': 'min'
        }).reset_index()
        flat_clusters = pd.concat([flat_clusters, cc])


    cluster_dict = quality_score(flat_clusters, running_difference)
    print(json.dumps(cluster_dict, indent=4))

    best_cluster_QS = -np.inf
    best_detection_dict = {}

    if cluster_dict:
        for d in cluster_dict:
            if cluster_dict[d]['THETA'] > best_cluster_QS:
                best_cluster_QS = cluster_dict[d]['THETA']
                best_detection_dict = cluster_dict[d]
    else:
        raise ValueError("Cluster dictionary empty")

    return best_detection_dict


# ─────────────────────────────────────────────
#  Entry point
# ─────────────────────────────────────────────

if __name__ == '__main__':

    data_set_testing = input('Want to test pre-downloaded image sets? (y/n): ')

    if data_set_testing.lower() == 'n':
        # ── Normal mode: optionally download, then run once ──
        decision = input('Want to download images? (y/n): ')
        if decision.lower() == 'y':
            folder = './data_processed/lasco/c3/'
            for filename in os.listdir(folder):
                file_path = os.path.join(folder, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('Failed to delete %s. Reason: %s' % (file_path, e))

            init_date = input('Insert real onset event date (YYYY-MM-DD HH:MM:SS): ')
            start_date, stop_date = get_date_range(init_date)

            download_metadata(start_date, folder)
            lasco = Lasco(start_date, stop_date)
            lasco.extract_data("c3")

        
        try:
            run_pipeline('./data_processed/lasco/c3')
        except e:
            print(e)

    else:
        # ── Test mode: loop over all pre-downloaded directories ──
        #
        # List every directory path you want to process here.
        # Each directory must contain PNGs + images_metadata.json.
        #
        test_dirs = []

        for dr in os.listdir('./data_processed/testing_images_c3/'):
            test_dirs.append('./data_processed/testing_images_c3/' + dr)

        print(test_dirs)

        all_results = {}
        for dir_path in test_dirs:
            if not os.path.isdir(dir_path):
                print(f"[SKIP] Directory not found: {dir_path}")
                continue
            result = run_pipeline(dir_path)
            all_results[dir_path] = result
            

        print("\n\n=== SUMMARY OF ALL RUNS ===")
        print(json.dumps(all_results, indent=4))

        real_speeds = {
            '20150621': 1366,
            '20040725': 1366,
            '20230421': 1284,
            '20231128': 741,
            '20240509': 1280,
            '20240808': 789,
            '20251216': 579,
            '20241009': 1435,
            '20250413': 529
        }

        s = 0
        l = 0
        for res in all_results.values():
            speed = res['velocity_km_s']
            date = res['onset_date']

            if date in real_speeds:
                s += abs(real_speeds[date] - speed)
                l += 1

                print(f'date {date}, real speed: {real_speeds[date]}, predicted speed: {speed}')

        print('MAE: ', s/l)