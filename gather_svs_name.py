import os
import csv
import glob

def get_files_with_extensions(root_dir, extensions):
    count = 1
    file_list = []
    for root, dirs, files in os.walk(root_dir):
        if 'wsi' in dirs:
            wsi_dir = os.path.join(root, 'wsi')
            os.chdir(wsi_dir)
            for ext in extensions:
                for file in glob.glob(f"*.{ext}"):
                    print(count)
                    count += 1
                    file_path = os.path.join(wsi_dir, file)
                    file_list.append(file_path)
    return file_list

def save_to_tsv(file_list, output_file):
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\n')
        writer.writerow(file_list)
    print(f"File list saved to '{output_file}'.")

# Example usage
root_directory = '/mnt/datas1/'  # Replace with the desired root directory
output_filename = '/mnt/datas2/histotil/svs_files.tsv'
extensions = ['svs', 'tif', 'tiff']  # Add or remove extensions as needed

files_with_extensions = get_files_with_extensions(root_directory, extensions)
save_to_tsv(files_with_extensions, output_filename)