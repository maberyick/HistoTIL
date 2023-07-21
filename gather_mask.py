import os
import shutil
import csv
import argparse

def find_and_copy_files(source_folder, destination_folder, extensions):
    copied_files = []
    for root, dirs, files in os.walk(source_folder):
        for ext in extensions:
            for file in files:
                if file.endswith(f'.{ext}_mask_use.png'):
                    source_path = os.path.join(root, file)
                    new_file_name = file.replace(f'.{ext}_mask_use.png', 'png')
                    destination_path = os.path.join(destination_folder, new_file_name)
                    shutil.copy2(source_path, destination_path)
                    copied_files.append(new_file_name)
                    print(f"Copied {file} to {destination_path}")

    return copied_files

def save_to_csv(file_names, csv_file):
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['File Name'])
        for file_name in file_names:
            writer.writerow([file_name])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find and copy files and save to CSV.")
    parser.add_argument("--source_folder", type=str, help="Source folder path")
    parser.add_argument("--destination_folder", type=str, help="Destination folder path")
    parser.add_argument("--csv_file", type=str, help="CSV file path")
    parser.add_argument("--extensions", nargs='+', default=['svs', 'tif', 'tiff'],
                        help="List of file extensions to check for (space-separated). Default: svs tif tiff")
    args = parser.parse_args()

    source_folder = args.source_folder
    destination_folder = args.destination_folder
    csv_file = args.csv_file
    extensions = args.extensions

    copied_files = find_and_copy_files(source_folder, destination_folder, extensions)
    save_to_csv(copied_files, csv_file)