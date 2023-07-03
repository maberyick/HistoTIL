import os
import shutil
import csv

def find_and_copy_files(source_folder, destination_folder):
    copied_files = []
    for root, dirs, files in os.walk(source_folder):
        for file in files:
            if file.endswith('.svs_mask_use.png'):
                source_path = os.path.join(root, file)
                new_file_name = file.replace('.svs_mask_use.png', '.png')
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

# Example usage
source_folder = r'E:\Tools\histoqc\histoqc_output_20230614-141029'
destination_folder = r'E:\wsi_data\tissue_mask'
csv_file = r'E:\wsi_data\tissue_mask_filename_list.csv'

copied_files = find_and_copy_files(source_folder, destination_folder)
save_to_csv(copied_files, csv_file)