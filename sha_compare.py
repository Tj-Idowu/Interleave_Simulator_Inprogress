#!/usr/bin/env python3

import hashlib
import sys
import os

# Nextflow variables (captured via template)
target_file = "${all_files}"
# We use .name if we want to match the filename string in the manifest
filename_to_match = "${all_files.name}"
manifest_path = "${checksum_sha}"

def verify():
    expected_hash = None
    
    # 1. Parse the standard .sha256 manifest
    # Format is: [hash][space][space][relative/path]
    try:
        with open(manifest_path, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.strip():
                    continue
                
                parts = line.split(None, 1) # Split on whitespace
                if len(parts) < 2:
                    continue
                
                current_hash = parts[0].strip().lower()
                current_path = parts[1].strip()
                
                # Match against the filename (or relative path)
                # We check if it ends with filename_to_match to handle relative paths
                if current_path == filename_to_match or current_path.endswith('/' + filename_to_match):
                    expected_hash = current_hash
                    break
    except FileNotFoundError:
        print(f"ERROR: Manifest file {manifest_path} not found.")
        sys.exit(1)

    if not expected_hash:
        print(f"ERROR: {filename_to_match} not found in manifest.")
        sys.exit(1)

    # 2. Calculate the SHA256 hash of the actual file
    sha256_hash = hashlib.sha256()
    try:
        with open(target_file, "rb") as f:
            # 64KB chunks for speed and memory efficiency
            for byte_block in iter(lambda: f.read(65536), b""):
                sha256_hash.update(byte_block)
    except FileNotFoundError:
        print(f"ERROR: File {target_file} not found on disk.")
        sys.exit(1)
    
    calculated_hash = sha256_hash.hexdigest().lower()

    # 3. Compare
    if calculated_hash == expected_hash:
        print(f"SUCCESS: {filename_to_match} integrity verified.")
        sys.exit(0)
    else:
        print(f"FAIL: Integrity Mismatch for {filename_to_match}!")
        print(f"Expected: {expected_hash}")
        print(f"Found:    {calculated_hash}")
        sys.exit(1)

if __name__ == "__main__":
    verify()