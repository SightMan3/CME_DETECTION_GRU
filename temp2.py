from datetime import datetime, timedelta
import os
import requests
import json
import shutil
from extractDataSrc import Lasco


def download_metadata(folder, time_str):
    original_datetime = datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S")
    date_str = original_datetime.strftime("%Y-%m-%dT%H:%M:%S.000Z")
    url = f"https://api.helioviewer.org/v2/getClosestImage/?date={date_str}&sourceId=5"

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad status codes
        
        data = response.json()
        
        # Save the JSON data to the specified path
        with open(f"{folder}/images_metadata" + ".json", 'w') as f:
            json.dump(data, f, indent=4)
    

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from Helioviewer API: {e}")
    except IOError as e:
        print(f"Error saving file to {folder}/images_metadata.json: {e}")

def get_date_range(date_str):
    dt = datetime.strptime(date_str, '%Y-%m-%d %H:%M:%S')
    
    start = (dt - timedelta(hours=2)).replace(minute=0, second=0)
    end = dt + timedelta(hours=12)
    
    return (
        start.strftime('%Y-%m-%d %H:%M:%S'),
        end.strftime('%Y-%m-%d %H:%M:%S')
    )

def download_images(folder, init_date):
    print('Downloading coronograph images -----------------')
        
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

    start_date, stop_date = get_date_range(init_date)

    download_metadata(folder, start_date)
    lasco = Lasco.Lasco(start_date, stop_date)
    lasco.extract_data("c3", custom_target_dir=folder)


if __name__ == '__main__':
    folder = './data_processed/testing_images_c3'
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

    dates = [
        '2015-06-21 02:35:00',
        '2000-04-04 16:43:00',
        '2004-07-25 14:54:00',
        '2023-04-21 18:12:00',
        '2023-11-28 20:12:00',
        '2024-05-09 09:24:00',
        '2024-08-08 19:48:00',
        '2025-12-16 09:34:00',
        '2024-10-09 02:12:00',
        '2025-04-13 06:36:00',
        '2025-01-04 05:48:00',
        '2025-01-04 18:36:00'
    ]

    folder = './data_processed/testing_images_c3'
    for i, date in enumerate(dates):
        os.mkdir(folder + f'/event{i}')

        download_images(folder + f'/event{i}', date)
    print(os.listdir(folder))
