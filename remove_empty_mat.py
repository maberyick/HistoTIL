import os
import argparse

def delete_small_mat_files(directory, size_limit_kb):
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".mat"):
                filepath = os.path.join(root, file)
                if os.path.getsize(filepath) < size_limit_kb * 1024:
                    print(f"Deleting: {filepath}")
                    os.remove(filepath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Delete small .mat files")
    parser.add_argument("path", type=str, help="Base directory path")
    parser.add_argument("--size-limit", type=int, default=1, help="Size limit in KB (default: 1)")

    args = parser.parse_args()
    base_directory = args.path
    size_limit_kb = args.size_limit
    target_directory = os.path.join(base_directory, "dataset_output")

    if os.path.exists(target_directory):
        delete_small_mat_files(target_directory, size_limit_kb)
        print("Deletion of small .mat files completed.")
    else:
        print("Target directory not found.")
