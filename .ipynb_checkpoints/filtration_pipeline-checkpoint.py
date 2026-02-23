import os, shutil
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
import cv2



def get_coronal_holes(start, end):
    '''
    ## Downloading data ##
    '''
    print('Downloading EIT data -----------------')
    start_eit_time = Time(start, scale='utc')
    end_eit_time = Time(end, scale='utc')
    
    time_range = a.Time(start, end)
    result = sunpy.net.Fido.search(
        time_range,  # start and end time
        a.Instrument.aia,  # SDO/AIA instrument
        a.Wavelength(193*u.angstrom),  # EUV 193 Å channel
        a.Sample(5*u.minute) 
    )
    
    # Download the files
    downloaded_files = sunpy.net.Fido.fetch(result)
    
    print(downloaded_files)


    '''
    ## Image augmentation ##
    '''

    SDO_map = sunpy.map.Map(downloaded_files[1])

    SDO_img = SDO_map.data
    
    # resized_eit = ndi.zoom(eit_data, IMG_SIZE / max(eit_data.shape), order=1)
    # normalized_eit = resized_eit / np.max(resized_eit)
    # eit_map_array_image = np.array(normalized_eit)
    
    
    '''
    Normalize brigthness
    '''
    vmin = 0
    vmax = 300
    bright_SDO = np.clip(SDO_img, vmin, vmax)  # keep values within [vmin, vmax]
    bright_SDO = ((bright_SDO - vmin) / (vmax - vmin) * 255).astype(np.uint8)
    
    less_bright_SDO = np.clip(SDO_img, vmin, vmax + 200)
    less_bright_SDO = ((less_bright_SDO - vmin) / ((vmax + 200) - vmin) * 255).astype(np.uint8)
    
    '''
    Crop image
    '''
    
    h, w = bright_SDO.shape
    cx, cy = w // 2, h // 2   # center (SOHO images are centered already)
    radius = min(cx, cy) - 400   # approximate disk radius
    
    # make a circular mask
    Y, X = np.ogrid[:h, :w]
    dist = np.sqrt((X - cx)**2 + (Y - cy)**2)
    mask = dist <= radius
    
    # apply mask
    cropped_bright_SDO = np.full_like(bright_SDO, 255)  # fill with white background
    cropped_bright_SDO[mask] = bright_SDO[mask]

    cropped_less_bright_SDO = np.full_like(less_bright_SDO, 255)  # fill with white background
    cropped_less_bright_SDO[mask] = less_bright_SDO[mask]
    
    '''
    Resize image
    
    CRB - cropped, resized, bright
    '''
    CRB_SDO = cv2.resize(cropped_bright_SDO, (256, 256), interpolation=cv2.INTER_AREA)
    
    
    # normalised, cropped, resized, brighthened SDO, NCRB - SDO
    pmin, pmax = np.percentile(CRB_SDO, (1, 99))
    NCRB_SDO = np.clip((CRB_SDO - pmin) / (pmax - pmin), 0, 1)
    
    NCRB_SDO = NCRB_SDO.reshape((1, 256, 256, 1))

    
    '''
    Load model and predict CHs
    '''
    from model_scss_net import scss_net
    from metrics import dice_np, iou_np, dice, iou
    from utils import plot_imgs, plot_metrics
    
    model = scss_net( 
        (256, 256, 1),
        filters=32,       
        layers=4,
        batch_norm=True,
        drop_prob=0.5)
    
    # Compile model
    model.compile(
        optimizer="adam",
        loss="binary_crossentropy",  
        metrics=[iou, dice])

    model.load_weights("./model_ch_region.h5")
    y_pred = model.predict(NCRB_SDO)

    pred_mask_resized = cv2.resize(y_pred[0], (4096, 4096), interpolation=cv2.INTER_NEAREST)

    plt.imshow(less_bright_SDO, cmap='grey')
    plt.imshow(pred_mask_resized, alpha=0.5)

    

def get_timeseries_data(start, end):
    '''
    Download Bx, By, Bz data and prepare
    '''
    print('Downloading B data -----------------')
    insitu = InSitu(start, end)
    insitu.extract_data("B", verbose=True)

    B = pd.read_csv('./data_processed/in_situ/B_data.csv')
    B.index = pd.to_datetime(B['Time'])
    B.drop('Time', axis=1, inplace=True)
    B['B'] = np.sqrt(B['B_x']**2 + B['B_y']**2 + B['B_z']**2) 

    B = B.resample('3min').mean()
    
    '''
    Download SYM_H data and prepare
    '''
    print('Downloading SYM_H -----------------')
    server = 'https://cdaweb.gsfc.nasa.gov/hapi'
    dataset = 'OMNI_HRO_1MIN'
    parameter = 'SYM_H'
    
    hapi_data = get_hapi_data(server, dataset, parameter, start, end)

    SYM_data = []
    for record in hapi_data:
        row = []
        row.append(record[0])
        row.append(round(record[1], 3))
        SYM_data.append(row)
    
    SYM_H = pd.DataFrame(data=SYM_data, columns=["Time", 'SYM_H'])
    SYM_H['Time'] = SYM_H['Time'].str.decode('utf-8')
    SYM_H.index = pd.to_datetime(SYM_H['Time'], format="%Y-%m-%dT%H:%M:%S.000Z")
    SYM_H.drop('Time', axis=1, inplace=True)
    SYM_H = SYM_H.resample('3min').mean()


    '''
    Get plasma temparature and beta data
    '''
    hapiT = get_hapi_data(server, 'OMNI_HRO2_1MIN', 'T,Beta', start, end)
    T = []
    for record in hapiT:
        row = []
        row.append(record[0])
        row.append(round(record[1], 3))
        row.append(round(record[2], 3))
        T.append(row)
    
    T = pd.DataFrame(data=T, columns=["Time", 'T', 'Beta'])
    T['Time'] = T['Time'].str.decode('utf-8')
    T.index = pd.to_datetime(T['Time'], format="%Y-%m-%dT%H:%M:%S.000Z")
    T.drop('Time', axis=1, inplace=True)
    T['T'] = T['T'].replace([999.9, 9999.9, 9999999.0], np.nan)
    T.loc[np.isclose(T["Beta"], 999.99), "Beta"] = np.nan
    T = T.resample('3min').mean()


    '''
    Calculate beta data
    '''
    k_B = 1.3806e-23
    mu_0 = 4 * np.pi * 1e-7
    
    
    '''
    Get the time range for imaging
    '''
    ts = SYM_H.idxmin().iloc[0]
    after = ts + timedelta(hours=3)
    before = ts - timedelta(hours=72)
    after, before

    return SYM_H, B, T, after, before
    

def get_Corongraph_InSitu(start, end):

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
    
    
    lasco = Lasco(start, end)
    lasco.extract_data("c3")

    # img_name_change = None
    
    # image_folder = './data_processed/lasco/c2/'
    # images = sorted(os.listdir(image_folder))
    # prev_img = None
    # for i, filename in enumerate(images):
    #     img_path = os.path.join(image_folder, filename)
    #     img = cv2.imread(img_path)

    #     if prev_img is not None:
    #         percent, _ = compute_brightness_change(img, prev_img)
    #         # print(f"percent of pixels changed for current image {filename} and before image : {percent}")

    #         if percent > 10:
    #             img_name_change = filename
    #             break
        
    #     prev_img = img

    # #20031028_1514_lascoc2_1024
    # if img_name_change is not None:
    #     img_split = img_name_change.split('_')
    #     formatted_YMD = f"{img_split[0][:4]}-{img_split[0][4:6]}-{img_split[0][6:]}"
    #     formatted_time = f"{img_split[1][:2]}:{img_split[1][2:4]}:00"
    #     full = f"{formatted_YMD} {formatted_time}"
    #     print(f"date of CME eruption {full}")

    #     return full
    # return None


Sun_struct_start = input('Sun structure date start: ')
Sun_struct_end   = str(datetime.strptime(Sun_struct_start, "%Y-%m-%d %H:%M:%S") + timedelta(hours=30))

Event_start = Sun_struct_end
Event_end   = input('Sun event end date: ')


SYM_H, B, T, after, before = get_timeseries_data(Event_start, Event_end)


try:
    get_Corongraph_InSitu(Sun_struct_start, Sun_struct_end)
except:
    print('Cant download LASCO c3 images')
    
try:
    get_coronal_holes(Sun_struct_start, 
                      str(datetime.strptime(Sun_struct_start, "%Y-%m-%d %H:%M:%S") + timedelta(minutes=20))
                     )
except:
    print('Cant download EIT images')

fig, ax = plt.subplots(3)

ax[0].plot(B['B_x'], label='$\\mathbf{B}_x$')
ax[0].plot(B['B_y'], label='$\\mathbf{B}_y$')
ax[0].plot(B['B_z'], label='$\\mathbf{B}_z$')
ax[0].plot(SYM_H, label='$SYM_H$')
ax[1].plot(T['T'], label='$\\mathbf{T}$')
ax[2].plot(T['Beta'], label='$\\beta$')
ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.show()


