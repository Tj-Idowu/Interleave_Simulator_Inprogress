##############################################
# .csv output file

import tarfile
import hashlib
import csv
import os

def calculate_tar_hashes(tar_path, output_csv):
    """
    Calculates SHA256 hashes for files inside a .tar.gz without extracting them.
    """
    # SHA256 block size (64kb) for efficient reading
    BLOCK_SIZE = 65536 

    print(f"Opening {tar_path}...")

    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            with open(output_csv, mode='w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                # Write header
                writer.writerow(['Filename', 'SHA256'])

                # Iterate through each item in the tarball
                for member in tar.getmembers():
                    # Only process regular files (skip directories and links)
                    if member.isfile():
                        sha256_hash = hashlib.sha256()
                        f = tar.extractfile(member)
                        
                        if f:
                            # Read the file in chunks to save memory
                            while True:
                                data = f.read(BLOCK_SIZE)
                                if not data:
                                    break
                                sha256_hash.update(data)
                            
                            # Record the result
                            writer.writerow([member.name, sha256_hash.hexdigest()])
                            print(f"Processed: {member.name}")

        print(f"\nSuccess! Hash list saved to: {output_csv}")

    except FileNotFoundError:
        print("Error: The specified tar.gz file was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # --- CONFIGURATION ---
    input_tarball = "your_data_file.tar.gz"  # input file path
    output_file = "file_hashes.csv"	# output file name
    # ---------------------

    calculate_tar_hashes(input_tarball, output_file)


##########################################################################################################################
# .sha output file


import tarfile
import hashlib
import os

def create_sha256_manifest(tar_path, output_manifest):
    """
    Calculates SHA256 hashes for files in a .tar.gz and writes a standard .sha256 manifest file.
    """
    BLOCK_SIZE = 65536 

    print(f"Reading archive: {tar_path}")

    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            with open(output_manifest, mode='w', encoding='utf-8') as f_out:
                
                # Iterate through members to get relative paths
                for member in tar.getmembers():
                    # Only process regular files
                    if member.isfile():
                        sha256_hash = hashlib.sha256()
                        
                        # extractfile gives us a stream of the file content
                        with tar.extractfile(member) as f_in:
                            while True:
                                data = f_in.read(BLOCK_SIZE)
                                if not data:
                                    break
                                sha256_hash.update(data)
                        
                        # Format: hash  relative/path/to/file
                        # Note: two spaces between hash and path is the standard
                        manifest_line = f"{sha256_hash.hexdigest()}  {member.name}\n"
                        f_out.write(manifest_line)
                        print(f"Hashed: {member.name}")

        print(f"\nSuccess! Standard manifest created: {output_manifest}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # --- CONFIGURATION ---
    input_tarball = "datasets1_10.tar.gz"
    output_manifest = "dataset1_10checksums.sha256"
    # ---------------------

    create_sha256_manifest(input_tarball, output_manifest)
