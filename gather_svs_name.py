import os
import csv
import glob
import argparse

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get files with extensions and save to TSV.")
    parser.add_argument("--root_directory", type=str, help="Root directory path")
    parser.add_argument("--output_filename", type=str, help="Output TSV file path")
    parser.add_argument("--extensions", nargs='+', default=['svs', 'tif', 'tiff'],
                        help="List of file extensions to search for (space-separated). Default: svs tif tiff")
    args = parser.parse_args()

    root_directory = args.root_directory
    output_filename = args.output_filename
    extensions = args.extensions

    files_with_extensions = get_files_with_extensions(root_directory, extensions)
    save_to_tsv(files_with_extensions, output_filename)
