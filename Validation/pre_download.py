import sys
from pathlib import Path

PROJECT_ROOT = Path.cwd().parent
sys.path.insert(0, str(PROJECT_ROOT))

from datetime import datetime, timedelta


import os
import json
import requests
import pandas as pd
from tqdm import tqdm
from extractDataSrc.Lasco import Lasco

CSV      = 'sample_w_CPA.csv'
OUT_ROOT = '/home/lukas/projects/CME_DETECTION_GRU/Validation/testing_images_c3'


def get_date_range(date_str):
    dt = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S')
    start = (dt - timedelta(hours=2)).replace(minute=0, second=0)
    end   = dt + timedelta(hours=12)
    return (
        start.strftime('%Y-%m-%d %H:%M:%S'),
        end.strftime('%Y-%m-%d %H:%M:%S')
    )


def download_metadata(time_str, save_dir):
    original_datetime = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
    date_str = original_datetime.strftime("%Y-%m-%dT%H:%M:%S.000Z")
    url = f"https://api.helioviewer.org/v2/getClosestImage/?date={date_str}&sourceId=5"

    response = requests.get(url, timeout=30)
    response.raise_for_status()
    with open(os.path.join(save_dir, "images_metadata.json"), 'w') as f:
        json.dump(response.json(), f, indent=4)


def download_image(init_date, folder):
    os.makedirs(folder, exist_ok=True)
    start_date, stop_date = get_date_range(init_date)

    download_metadata(start_date, folder)
    lasco = Lasco(start_date, stop_date)
    lasco.extract_data(custom_target_dir=folder, detector="c3", skip_images=0, verbose=False)


if __name__ == '__main__':
    cols = ['index', 'date', 'time', 'central_pa', 'mpa', 'width', 'linear_speed', 'mass', 'mass_flagged']
    events = pd.read_csv(CSV, names=cols, header=0)

    for _, row in tqdm(events.iterrows(), total=len(events), unit='event'):
        init_date = f"{row['date'].replace('/', '-')} {row['time']}"
        dt = datetime.strptime(init_date, '%Y-%m-%d %H:%M:%S')
        folder = os.path.join(OUT_ROOT, f"{int(row['index'])}_{dt.strftime('%Y%m%d_%H%M%S')}")

        # already done -> skip, so you can re-run after a crash
        if os.path.exists(os.path.join(folder, 'event.json')):
            continue

        try:
            download_image(init_date, folder)
            with open(os.path.join(folder, 'event.json'), 'w') as f:
                json.dump({'date': init_date, 'linear_speed': row['linear_speed'],
                           'width': row['width'], 'mass': row['mass']}, f, indent=4)
        except Exception as e:
            tqdm.write(f"[FAIL] {init_date}: {e}")