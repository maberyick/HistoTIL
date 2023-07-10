import os
import csv
import glob

def get_svs_files(root_dir):
    count = 1
    file_list = []
    for root, dirs, files in os.walk(root_dir):
        if 'wsi' in dirs:
            wsi_dir = os.path.join(root, 'wsi')
            os.chdir(wsi_dir)
            for file in glob.glob("*.svs"):
                if file.endswith('.svs'):
                    print(count)
                    count = count+1
                    file_path = os.path.join(wsi_dir, file)
                    file_list.append((file_path))
    return file_list

def save_to_tsv(file_list, output_file):
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\n')
        writer.writerow(file_list)
    print(f"File list saved to '{output_file}'.")

# Example usage
root_directory = '/mnt/datas1/'  # Replace with the desired root directory
output_filename = '/mnt/datas2/histotil/svs_files.tsv'

svs_files = get_svs_files(root_directory)
save_to_tsv(svs_files, output_filename)