import os
import argparse

def is_empty_folder(folder_path):
    return not any(os.scandir(folder_path))

def delete_empty_subfolders(root_path):
    for folder_name in os.listdir(root_path):
        folder_path = os.path.join(root_path, folder_name, 'TIL_features')
        folder_path_to_delete = os.path.join(root_path, folder_name)
        if os.path.isdir(folder_path) and is_empty_folder(folder_path):
            print(f"Deleting empty subfolder: {folder_path}")
            os.rmdir(folder_path)
            os.rmdir(folder_path_to_delete)

def main():
    parser = argparse.ArgumentParser(description="Delete empty subfolders in a given path.")
    parser.add_argument("path", metavar="PATH", type=str, help="The path to the parent folder.")
    args = parser.parse_args()

    target_path = args.path

    if not os.path.exists(target_path):
        print(f"The provided path '{target_path}' does not exist.")
    else:
        delete_empty_subfolders(target_path)

if __name__ == "__main__":
    main()
