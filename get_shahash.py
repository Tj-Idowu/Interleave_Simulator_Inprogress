#!/usr/bin/env python3

import tarfile
import hashlib
import os
import glob


def create_sha256_manifest(tar_files, output_manifest):

    BLOCK_SIZE = 65536  # Read files in 64 KB blocks

    try:
        # open a manifest for all .tar.gz files in the input folder
        with open(output_manifest, mode="w", encoding="utf-8") as f_out:

            # process each file
            for tar_path in tar_files:
                print(f"\nReading archive: {tar_path}")
                try:
                    with tarfile.open(tar_path, "r:gz") as tar:
                        # process each file inside the archive
                        for member in tar:
                            if member.isfile():
                                sha256_hash = hashlib.sha256()
                                # open the archived file without extracting it
                                f_in = tar.extractfile(member)
                                if f_in is None:
                                    raise OSError(
                                        f"Could not read archived file: "
                                        f"{member.name}"
                                    )
                                # read and hash the file in blocks
                                with f_in:
                                    while True:
                                        data = f_in.read(BLOCK_SIZE)
                                        if not data:
                                            break
                                        sha256_hash.update(data)
                                # include archive name to the manifest
                                archive_name = os.path.basename(tar_path)
                                manifest_line = (
                                    f"{sha256_hash.hexdigest()}  "
                                    f"{archive_name}::{member.name}\n"
                                )
                                f_out.write(manifest_line)
                                print(f"Hashed: {archive_name}::{member.name}")
                except (tarfile.TarError, OSError) as e:
                    print(f"Error processing archive {tar_path}: {e}")
                except Exception as e:
                    print(
                        f"Unexpected error processing archive "
                        f"{tar_path}: {e}"
                    )
        print(f"\nAll available checksums written to: {output_manifest}")
    except OSError as e:
        print(f"Error creating manifest {output_manifest}: {e}")

if __name__ == "__main__":

    # input folder
    input_folder = "datasets1_10"
    # find all .tar.gz files in the folder
    tar_files = sorted(
        glob.glob(os.path.join(input_folder, "*.tar.gz"))
    )
    # output manifest file name
    output_manifest = os.path.join(
        input_folder,
        "dataset_checksums.sha256"
    )
    if len(tar_files) == 0:
        print(f"Error: No .tar.gz files found in {input_folder}")
    else:
        print(f"Found {len(tar_files)} archive(s) in {input_folder}")
        # calculate all checksums and write them to .sha256 manifest file
        create_sha256_manifest(tar_files, output_manifest)
