import os, shutil

import requests
from extractDataSrc.Lasco import Lasco
from extractDataSrc.InSitu import InSitu
from extractDataSrc.Eit195 import Eit195
from datetime import datetime, timedelta

import sunpy.timeseries as ts
from sunpy.net import Fido, attrs as a
from extractDataSrc.Hapi import get_hapi_data

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sunpy
from sunpy.net import attrs as a
from astropy.time import Time
import astropy.units as u
import sunpy.map
import scipy.ndimage as ndi
import albumentations
import cv2 as cv
import scipy

from PIL import Image
import json


from skimage.measure import label, regionprops
import math
from matplotlib.patches import Ellipse

from skimage import exposure
from skimage.feature import canny

from skimage.transform import hough_line, hough_line_peaks
from skimage.filters import meijering
from skimage.morphology import skeletonize
from skimage.morphology import disk, binary_opening


decision = input('Want to dowload images? (y/n) ')
if decision.lower() == 'y':
    
    """
        Loading images and metadata
    """
    def download_metadata(time_str):
        original_datetime = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
        date_str = original_datetime.strftime("%Y-%m-%dT%H:%M:%S.000Z")
        url = f"https://api.helioviewer.org/v2/getClosestImage/?date={date_str}&sourceId=5"

        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an error for bad status codes
            
            data = response.json()
            
            # Save the JSON data to the specified path
            with open("./data_processed/lasco/c3/images_metadata" + ".json", 'w') as f:
                json.dump(data, f, indent=4)
        

        except requests.exceptions.RequestException as e:
            print(f"Error fetching data from Helioviewer API: {e}")
        except IOError as e:
            print(f"Error saving file to ./data_processed/lasco/c3/images_metadata.json: {e}")


    print('Downloading coronograph images -----------------')
        
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

    start_date = '2015-06-21 00:00:00'
    stop_date = '2015-06-21 10:00:00'

    download_metadata(start_date)
    lasco = Lasco(start_date, stop_date)
    lasco.extract_data("c3")


images = []

for im in os.listdir('./data_processed/lasco/c3'):
    if (im.endswith(".png")):
        img = Image.open('./data_processed/lasco/c3/' + im)
        images.append(np.array(img.convert('L')))


png_file = f'data_processed/lasco/c3/{os.listdir("data_processed/lasco/c3")[0]}'
img = Image.open(png_file)
image_data = np.array(img.convert('L'))

json_file = 'data_processed/lasco/c3/images_metadata.json'

with open(json_file, 'r') as f:
    metadata = json.load(f)



center_x = 512
center_y = 512

arcsec_per_pixel = metadata['scale']
ARCSEC_TO_KM = 725.0    # Approximate conversion at 1 AU
V_MIN_KM_S = 100.0      # Minimum plausible CME speed
V_MAX_KM_S = 2000.0     # Maximum plausible CME speed
rsun_pixels = metadata['rsun']

height, width = image_data.shape

y_indices, x_indices = np.indices((height, width))

x_rel_pixels = x_indices - center_x
y_rel_pixels = y_indices - center_y

x_arcsec = x_rel_pixels * arcsec_per_pixel
y_arcsec = y_rel_pixels * arcsec_per_pixel

r_arcsec = np.sqrt(x_arcsec**2 + y_arcsec**2)


"""
    Functions
"""

def image_masking(image_data):
    mask = (r_arcsec < 4.7*rsun_pixels*arcsec_per_pixel
       ) | (r_arcsec > 512*arcsec_per_pixel
       ) | ((theta_arcsec < (3/4)*np.pi) & (theta_arcsec > (3/4 - 7/100)*np.pi))

    image_data_masked = image_data.copy().astype(np.float32)
    image_data_masked[mask] = np.nan

    return image_data_masked



def prepare_image(image):
    height, width = image_data.shape

    y_indices, x_indices = np.indices((height, width))

    x_rel_pixels = x_indices - center_x
    y_rel_pixels = y_indices - center_y

    x_arcsec = x_rel_pixels * arcsec_per_pixel
    y_arcsec = y_rel_pixels * arcsec_per_pixel

    r_arcsec = np.sqrt(x_arcsec**2 + y_arcsec**2)

    """
    Mask out sun, occulter, edges
    """
    theta_arcsec = np.arctan2(y_arcsec, x_arcsec)
    mask = (r_arcsec < 4.7*rsun_pixels*arcsec_per_pixel
       ) | (r_arcsec > 512*arcsec_per_pixel
       ) | ((theta_arcsec < (3/4)*np.pi) & (theta_arcsec > (3/4 - 7/100)*np.pi))

    image_data_masked = image.copy().astype(np.float32)
    image_data_masked[mask] = np.nan

    return image_data_masked, r_arcsec, theta_arcsec



def unroll_to_polar(image_data, r_arcsec_array, theta_arcsec_array, 
                    num_radii=512, num_angles=360):
    """
    Converts a single Cartesian image (image_data) to a Polar (R, Theta) map.

    Args:
        image_data (np.ndarray): The 2D image data (after masking).
        r_arcsec_array (np.ndarray): The calculated radial distances (in arcsec).
        theta_arcsec_array (np.ndarray): The calculated angles (in radians).
        num_radii (int): The desired number of radial bins.
        num_angles (int): The desired number of angular bins (0-360 degrees).

    Returns:
        np.ndarray: The 2D polar map (num_radii x num_angles).
    """
    height, width = image_data.shape
    
    R_max = r_arcsec_array.max()
    R_new = np.linspace(r_arcsec_array.min(), R_max, num_radii)
    Theta_new = np.linspace(0, 2 * np.pi, num_angles, endpoint=False)
    R_grid, Theta_grid = np.meshgrid(R_new, Theta_new)
    
    occulter_center_x = 512.0
    occulter_center_y = 512.0
    
    indices_x = occulter_center_x + R_grid * np.cos(Theta_grid) / arcsec_per_pixel
    indices_y = occulter_center_y + R_grid * np.sin(Theta_grid) / arcsec_per_pixel

    coords = np.vstack([indices_y.ravel(), indices_x.ravel()])

    polar_map = scipy.ndimage.map_coordinates(
        input=image_data, 
        coordinates=coords, 
        order=1,  # Linear interpolation
        cval=np.nan # Use NaN for masked/out-of-bounds areas
    ).reshape(num_angles, num_radii)
    
    return polar_map.T



def create_jmap(diff_ims, angle):
    j_map_slices = []
    
    for polar_map in diff_ims[::-1]:
        radial_slice = polar_map[:, angle]
        j_map_slices.append(radial_slice)
        
    j_map = np.stack(j_map_slices, axis=0) 
    
    return j_map


def cme_blob_detection(j_map, verbose=True):
    j_map_clean = np.nan_to_num(j_map, nan=0.0, posinf=0.0, neginf=0.0)
    p1, p99 = np.percentile(j_map_clean, (1, 99))
    j_map_clipped = np.clip(j_map_clean, p1, p99)
    j_map_norm = (j_map_clipped - j_map_clipped.min()) / (j_map_clipped.max() - j_map_clipped.min())
    j_map_eq = exposure.equalize_adapthist(j_map_norm, clip_limit=0.03)

    j_map_blobs = meijering(j_map_eq, sigmas=range(1, 4), black_ridges=False)
    threshold = np.percentile(j_map_blobs, 95) 
    binary_j_map = j_map_blobs > threshold

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
        center_xy = (x0, y0)


        
        """
        Start leading edge calculation of slope
        """
        x_coords = props.coords[:, 1]
        y_coords = props.coords[:, 0]
        
        unique_radii = np.unique(x_coords)
        leading_edge_y_min = []
        leading_edge_y_max = []
        leading_edge_x = []

        min_y_edge = [np.min(y_coords[x_coords == r]) for r in unique_radii]
        max_y_edge = [np.max(y_coords[x_coords == r]) for r in unique_radii]

        # Fit both
        slope_min = np.polyfit(unique_radii, min_y_edge, 1)[0]
        slope_max = np.polyfit(unique_radii, max_y_edge, 1)[0]
    
        leading_edge_slope = min(abs(slope_min), abs(slope_max))


        """
        End leading edge slope calculation
        """

        # slope = math.tan(angle_rad)
        slope = leading_edge_slope
        if slope <= 0.02:
            continue

        R_arcsec_max = (1024 / 2) * arcsec_per_pixel
        R_arcsec_min = 3.7 * rsun_pixels * arcsec_per_pixel
        dr_arcsec_per_pixel = (R_arcsec_max - R_arcsec_min) / 512 
        dt_seconds = 12 * 60.0 

        velocity_arcsec_per_sec = abs(1.0 / slope) * (dr_arcsec_per_pixel / dt_seconds)
        velocity_km_s = velocity_arcsec_per_sec * ARCSEC_TO_KM
        
        t_onset_idx = abs(y0 - slope * x0)
        
        ellipse_patch = Ellipse(
            xy=center_xy,
            width=major_axis,
            height=minor_axis,
            angle=-math.degrees(angle_rad),
            edgecolor='red',   
            facecolor='none',  
            linewidth=2,
            linestyle='--'             
        )

        if V_MIN_KM_S < velocity_km_s < V_MAX_KM_S:
            cme_detection.append({
                'x0': x0,
                'y0': y0,
                'angle_deg': math.degrees(angle_rad),
                'major_axis': major_axis,
                'minor_axis': minor_axis,
                'slope': slope,
                'velocity_km_s': velocity_km_s,
                't_onset_idx': math.ceil(t_onset_idx)
            })

        if verbose and slope < 1:
            ax.add_patch(ellipse_patch)
            ax.plot(np.linspace(57, 100, 100), -slope*np.linspace(0, 100, 100) + y0)
            
            print(f"Blob at (x={x0:.0f}, y={y0:.0f}): Angle deg{math.degrees(angle_rad):.1f} deg, Length {major_axis:.1f}, angle rad {angle_rad:.1f}rad")
            print('slope: ', slope, 'velocity: ', velocity_km_s, 't_onset_index: ', t_onset_idx)
            print()


    return cme_detection


def all_angles_detections_blobs(images):
    rows = []
    for angle in range(0, 360, 5):
        j_map = create_jmap(images, angle)
        cme = cme_blob_detection(j_map, verbose=False)

        for r in cme:
            r['angle'] = angle
            rows.append(r)

    return pd.DataFrame(rows)



def fill_factor(cluster):
    c1_angles = cluster['angle'].sort_values().values # Ensure it's a numpy array
    
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



def quality_score(flat_clusters):
    cluster_dict = {}

    for cluster_id in flat_clusters['cme_cluster_id'].unique():

        cluster = flat_clusters[flat_clusters['cme_cluster_id'] == cluster_id]
        
        v = cluster['velocity_km_s']
        t = cluster['t_onset_idx']

        cv_t = t.std() / t.mean() if t.mean() != 0 else 1.0
        cv_v = v.std() / v.mean() if v.mean() != 0 else 1.0

        FF, AW = fill_factor(cluster)

        norm_AW = AW / 360.0
        norm_FF = FF / 100.0
        
        norm_cv_v = max(0.0, 1.0 - cv_v)
        norm_cv_t = max(0.0, 1.0 - cv_t)

        qs = (0.4 * norm_AW) + (0.3 * norm_FF) + (0.15 * norm_cv_v) + (0.15 * norm_cv_t)

        cluster_dict[cluster_id.item()] = {
            "zscore_velocity": cv_v.item(),
            "zscore_onset_time": cv_t.item(),
            "velocity_km_s": v.quantile(0.9).item(),
            "onset_time": t.mean().item(),
            "angular_width": AW.item(),
            'fill_factor': FF.item(),
            'THETA': qs.item()
        }

    return cluster_dict


"""
    Ellipse detection for images
"""

prepared_images = []
prepared_r_arcsecs = []
prepared_theta_arcsecs = []

for image in images:
    im, r, t = prepare_image(image)
    prepared_images.append(im)
    prepared_r_arcsecs.append(r)
    prepared_theta_arcsecs.append(t)



polar_maps = []

for image, r_arcsec, theta_arcsec in zip(prepared_images, prepared_r_arcsecs, prepared_theta_arcsecs):
    polar_maps.append(unroll_to_polar(image, r_arcsec, theta_arcsec))



running_difference = []

for i in range(1, len(polar_maps)):
    running_difference.append(
        np.clip(polar_maps[i] - polar_maps[i-1], 0, 255)
    )
    # running_difference.append(
    #     polar_maps[i] - polar_maps[i-1]
    # )





detections = all_angles_detections_blobs(running_difference)

detections = detections.sort_values(by='t_onset_idx').reset_index(drop=True)

detections['onset_diff'] = detections['t_onset_idx'].diff()
detections.fillna(0, inplace=True)
detections['cme_cluster_id'] = (detections['onset_diff'] > 1).cumsum()



counts = detections['cme_cluster_id'].value_counts()
valid_clusters = counts[counts > 10].index
detections = detections[detections['cme_cluster_id'].isin(valid_clusters)]

flat_clusters = pd.DataFrame()

for c in detections['cme_cluster_id'].unique():
    cc = detections[detections['cme_cluster_id'] == c].groupby('angle').agg({
            't_onset_idx': 'median', # Pandas automatically finds the median timestamp!
            'velocity_km_s': 'median',
            'cme_cluster_id': 'min'
        }).reset_index()

    flat_clusters = pd.concat([flat_clusters, cc])


cluster_dict = quality_score(flat_clusters)
print(json.dumps(cluster_dict, indent=4))


"""
Version 0.2
Added onset time calculation as t_onset_idx = abs(y0 - slope * x0) and added it to the cme_detection dictionary.
This should make onseet time approximation better
"""