#!/usr/bin/env python3

"""
Purpose
-------
This script takes a background metagenome file, a pathogen file, and its corresponding decoy file an shuffles them together.

Input
------
The inputs are the simulated metagenome, pathogen, and decoy reads as well as a number list for consistent randomisation.

Authorship
----------
Olateju Idowu
https://github.com/Tj-Idowu
"""

import os
import itertools
import gzip
import sys
import tempfile
import shutil
from Bio import SeqIO

# --- Configuration Paths ---
simulated_dir = "${params.bg_dir}"
patho_decoy_dir = "${params.patho_decoy_dir}"
number_list_file = "${params.numbers}"
final_output_dir = "${params.output_dir}"

# --- Function that interleaves the background metagenome with the pathogen and decoy files into a single file for each pathogen/decoy pair ---
def interleave_fastqs(file1, file2, file3, output_path):
    with gzip.open(output_path, 'wt') as out_handle:
        # use context managers for input files to make sure they close properly
        with gzip.open(file1, 'rt') as f1, \
             gzip.open(file2, 'rt') as f2, \
             gzip.open(file3, 'rt') as f3:
                
            records1 = SeqIO.parse(f1, 'fastq')
            records2 = SeqIO.parse(f2, 'fastq')
            records3 = SeqIO.parse(f3, 'fastq')
                
            # make sure all the files are processed regardless of their length
            for r1, r2, r3 in itertools.zip_longest(records1, records2, records3):
                if r1: SeqIO.write(r1, out_handle, 'fastq')
                if r2: SeqIO.write(r2, out_handle, 'fastq')
                if r3: SeqIO.write(r3, out_handle, 'fastq')

# --- Function to Read in The Numbers for The Sequences ---
def read_number_list(number_list_file):
    with open(number_list_file, 'r') as f:
        # read all lines, strip whitespace, and convert to integer
        # handle files where numbers might be on separate lines or space/comma separated
        numbers = [int(num) for line in f for num in line.strip().split() if num.isdigit()]
    return numbers

# --- Function to Determine Which Files to Interleave Together ---
def process_and_interleave():
    # create initial interleaved temporary files
    bg_map = {
        "Banthracis": "Simulated_WastewaterwoBacillus",
        "Cbotulinum": "Simulated_WastewaterwoClostridium",
        "Ecoli":      "Simulated_WastewaterwoEscherichia",
        "Ypestis":    "Simulated_WastewaterwoYersinia",
        "default":    "Simulated_WastewaterwoFrancisella"
    }

    all_files = os.listdir(patho_decoy_dir)
    patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq.gz")])
    decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq.gz")])

    # store the paths to the temporary files to pass to the next step
    interleaved_results = []

    for patho_file, decoy_file in zip(patho_files, decoy_files):
        bg_prefix = bg_map["default"]
        for key in bg_map:
            if key in patho_file:
                bg_prefix = bg_map[key]
                break

        # setup input paths (background, pathogens and decoys)
        patho_path_R1 = os.path.join(patho_decoy_dir, patho_file)
        patho_path_R2 = patho_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')
        decoy_path_R1 = os.path.join(patho_decoy_dir, decoy_file)
        decoy_path_R2 = decoy_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')
        sim_file_R1 = os.path.join(simulated_dir, f"{bg_prefix}_R1.fastq.gz")
        sim_file_R2 = os.path.join(simulated_dir, f"{bg_prefix}_R2.fastq.gz")

        output_prefix = "WW_" + patho_file.split('_', 1)[1].rsplit('_R1', 1)[0]

        tmp_R1 = tempfile.NamedTemporaryFile(prefix=f"{output_prefix}_", suffix="_R1.fastq.gz", delete=False)
        tmp_R2 = tempfile.NamedTemporaryFile(prefix=f"{output_prefix}_", suffix="_R2.fastq.gz", delete=False)
        tmp_R1.close()
        tmp_R2.close()

        interleave_fastqs(sim_file_R1, patho_path_R1, decoy_path_R1, tmp_R1.name)
        interleave_fastqs(sim_file_R2, patho_path_R2, decoy_path_R2, tmp_R2.name)

        interleaved_results.append({
            "sample_id": output_prefix,
            "r1_path": tmp_R1.name,
            "r2_path": tmp_R2.name
        })
    return interleaved_results

# --- Function to Shuffle The Interleaved Files ---
def shuffle_and_extract(data_list, number_list):
    # shuffles reads based on number_list into new temp files
    shuffled_results = []
    for entry in data_list:
        shuffled_R1 = tempfile.NamedTemporaryFile(prefix=f"shuffled_{entry['sample_id']}_", suffix="_R1.fastq.gz", delete=False)
        shuffled_R2 = tempfile.NamedTemporaryFile(prefix=f"shuffled_{entry['sample_id']}_", suffix="_R2.fastq.gz", delete=False)
        shuffled_R1.close()
        shuffled_R2.close()

        # index reads in memory
        read_map_R1, read_map_R2 = {}, {}
        with gzip.open(entry['r1_path'], 'rt') as f1, gzip.open(entry['r2_path'], 'rt') as f2:
            for i, (rec1, rec2) in enumerate(zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')), 1):
                read_map_R1[i], read_map_R2[i] = rec1, rec2

        # write in specified order
        with gzip.open(shuffled_R1.name, 'wt') as out1, gzip.open(shuffled_R2.name, 'wt') as out2:
            for num in number_list:
                if num in read_map_R1:
                    SeqIO.write(read_map_R1[num], out1, 'fastq')
                    SeqIO.write(read_map_R2[num], out2, 'fastq')

        shuffled_results.append({
            "sample_id": entry['sample_id'],
            "r1_path": shuffled_R1.name,
            "r2_path": shuffled_R2.name
        })
        # cleanup intermediate temporary files
        os.remove(entry['r1_path'])
        os.remove(entry['r2_path'])
        
    return shuffled_results

# --- Export shuffled files with subdirectory logic ---
def export_final_files(data_list, destination):
    # moves files from /tmp to permanent storage with folder-specific names
    for entry in data_list:
        sample_id = entry['sample_id']
        
        # determine the correct subdirectory based on the name
        sub_folder = ""
        if "_f1_" in sample_id.lower():
            sub_folder = "F1"
        elif "_f2_" in sample_id.lower():
            sub_folder = "F2"
        elif "_f3_" in sample_id.lower():
            sub_folder = "F3"
        
        # construct the full target path
        target_dir = os.path.join(destination, sub_folder)
        
        # create the subdirectory (and parent destination) if they don't exist
        os.makedirs(target_dir, exist_ok=True)
        
        final_R1 = os.path.join(target_dir, f"spiked_{sample_id}_R1.fastq.gz")
        final_R2 = os.path.join(target_dir, f"spiked_{sample_id}_R2.fastq.gz")
        
        # move the files from temp to the organized subfolders
        shutil.move(entry['r1_path'], final_R1)
        shutil.move(entry['r2_path'], final_R2)
        
        print(f"Exported: {final_R1}")

# --- Main Pipeline Execution ---
if __name__ == "__main__":
    # Interleave
    temp_interleaved = process_and_interleave()
    
    # Get Numbers
    order_numbers = read_number_list(number_list_file)
    
    # Shuffle
    temp_shuffled = shuffle_and_extract(temp_interleaved, order_numbers)
    
    # Export
    export_final_files(temp_shuffled, final_output_dir)