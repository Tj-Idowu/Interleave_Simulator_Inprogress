#!/usr/bin/env python3

# Packages used
import os
from Bio import SeqIO
import itertools
import gzip
import sys
from collections import defaultdict

# --- Configuration Paths ---
# Path to directories and files needed
simulated_dir = "/home/tj/Downloads/Benchmarking/Interleave_Simulator/test/bg_metagenome"
patho_decoy_dir = "/home/tj/Downloads/Benchmarking/Interleave_Simulator/test/patho_decoy"
output_dir = "/home/tj/Downloads/Benchmarking/Interleave_Simulator/test/interleaved"
number_list_file = "/home/tj/Downloads/Benchmarking/Interleave_Simulator/fastq_number_order.txt"

# Function to read in the numbers for the sequences
def read_number_list(number_list_file):
    """
    Read the input number list from the given text file.
    Return a list of integers.
    """
    try:
        with open(number_list_file, 'r') as f:
            # Read all lines, strip whitespace, and convert to integer
            # Handle files where numbers might be on separate lines or space/comma separated
            numbers = [int(num) for line in f for num in line.strip().split() if num.isdigit()]
        return numbers
    except FileNotFoundError:
        print(f"Error: {number_list_file} does not exist", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: cannot read {e}", file=sys.stderr)
        sys.exit(1)

# Function to number the sequences in the FASTQ files and reorder them based on the number list

def process_and_reorder_fastq(input_file, numbers_to_keep, temp_output_file):
    """
    1. Number all sequences in the input FASTQ file and stores them by number.
    2. Reorder the sequences based on the 'numbers_to_keep' list.
    3. Write the reordered sequences to a temporary output file. (merge 2 and 3)

    The FASTQ headers will be modified from:
    @<original_id> to: @<input_number> <original_id> 
    """
    print(f"Processing: {os.path.basename(input_file)}")
    
    # 1. Number all input sequences and store in a dictionary by their input number
    sequence_map = {}
    try:
        with gzip.open(input_file, 'rt') as in_handle:
            # Use enumerate to assign a 1-based input number
            for input_number, record in enumerate(SeqIO.parse(in_handle, 'fastq'), 1):
                # Store the *original* record, but modify the ID when writing out.
                sequence_map[input_number] = record
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        raise # Re-raise to be caught by the main logic
    except Exception as e:
        print(f"An unexpected error occurred during reading {input_file}: {e}", file=sys.stderr)
        raise

    # 2. Reorder sequences and write to temporary file
    # Use 'wt' to write a text file
    print(f"Reordering and writing to temporary file: {os.path.basename(temp_output_file)}")
    with open(temp_output_file, 'w') as out_handle:
        records_written = 0
        for target_number in numbers_to_keep:
            # Check if the target number is in the sequence map (i.e., if that sequence exists)
            if target_number in sequence_map:
                record = sequence_map[target_number]
                
                # Update the header to include the input number (the description is preserved)
                original_id = record.id
                original_description = record.description.split(' ', 1)[1] if ' ' in record.description else ''
                
                # New ID format: @<input_number> <original_description>
                record.id = str(target_number)
                record.description = f"{str(target_number)} {original_description}"
                
                SeqIO.write(record, out_handle, 'fastq')
                records_written += 1
            # else: skip the number if it's not in the file

    print(f"Finished processing. {len(sequence_map)} input sequences. {records_written} sequences written to the temporary file.")
    return temp_output_file


# --- Interleave Function (Modified to accept uncompressed temp files) ---

def interleave_fastqs(file1, file2, file3, output_file):
    """
    Interleaves three FASTQ files (file1, file2, file3) into a single gzipped output.
    All input files are now assumed to be uncompressed text files for simplicity,
    as the reorder function writes to uncompressed temp files.
    """
    print(f"Interleaving {os.path.basename(file1)}, {os.path.basename(file2)}, {os.path.basename(file3)} -> {os.path.basename(output_file)}")
    
    # Open the output file with gzip.open so that the output is in gz compressed format
    with gzip.open(output_file, 'wt') as out_handle:
        try:
            # Open input files as standard text files ('r')
            records1 = SeqIO.parse(open(file1, 'r'), 'fastq')
            records2 = SeqIO.parse(open(file2, 'r'), 'fastq')
            records3 = SeqIO.parse(open(file3, 'r'), 'fastq')
            
            # Make sure all files are processed regardless of length
            for r1, r2, r3 in itertools.zip_longest(records1, records2, records3):
                if r1:
                    SeqIO.write(r1, out_handle, 'fastq')
                if r2:
                    SeqIO.write(r2, out_handle, 'fastq')
                if r3:
                    SeqIO.write(r3, out_handle, 'fastq')
                    
        except FileNotFoundError as e:
            # This should not happen if the check below is done, but kept for robustness
            print(f"Error: One of the input files was not found during interleaving: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred during interleaving: {e}", file=sys.stderr)
            sys.exit(1)
            
    print(f"Successfully interleaved to {os.path.basename(output_file)}")


# --- Main Logic ---

def main():
    """Main function to orchestrate the process."""
    
    # Check that output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the master list of numbers ONCE
    numbers_to_keep = read_number_list(number_list_file)
    print(f"Loaded {len(numbers_to_keep)} target sequence numbers from {os.path.basename(number_list_file)}")

    # Get the pathogen and decoy files
    try:
        all_files = os.listdir(patho_decoy_dir)
    except FileNotFoundError:
        print(f"Error: Pathogen/Decoy directory not found at {patho_decoy_dir}", file=sys.stderr)
        sys.exit(1)
        
    patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq.gz")])
    decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq.gz")])

    if not patho_files or not decoy_files:
        print("Warning: No pathogen or decoy R1 files found. Exiting.", file=sys.stderr)
        return

    # Get the simulated wastewater files (assumed to be correct and pre-gzipped)
    sim_files = {
        "Bacillus": (os.path.join(simulated_dir, "Simulated_WastewaterwoBacillus_R1.fastq.gz"), 
                     os.path.join(simulated_dir, "Simulated_WastewaterwoBacillus_R2.fastq.gz")),
        "Clostridium": (os.path.join(simulated_dir, "Simulated_WastewaterwoClostridium_R1.fastq.gz"), 
                        os.path.join(simulated_dir, "Simulated_WastewaterwoClostridium_R2.fastq.gz")),
        "Escherichia": (os.path.join(simulated_dir, "Simulated_WastewaterwoEscherichia_R1.fastq.gz"), 
                        os.path.join(simulated_dir, "Simulated_WastewaterwoEscherichia_R2.fastq.gz")),
        "Yersinia": (os.path.join(simulated_dir, "Simulated_WastewaterwoYersinia_R1.fastq.gz"), 
                     os.path.join(simulated_dir, "Simulated_WastewaterwoYersinia_R2.fastq.gz")),
        "Francisella": (os.path.join(simulated_dir, "Simulated_WastewaterwoFrancisella_R1.fastq.gz"), 
                        os.path.join(simulated_dir, "Simulated_WastewaterwoFrancisella_R2.fastq.gz"))
    }


    # Process each pair of pathogen and decoy files
    for patho_file, decoy_file in zip(patho_files, decoy_files):
        patho_path_R1 = os.path.join(patho_decoy_dir, patho_file)
        decoy_path_R1 = os.path.join(patho_decoy_dir, decoy_file)
        
        # Determine the background file prefix
        if "Banthracis" in patho_file or "Banthracis" in decoy_file:
            sim_key = "Bacillus"
        elif "Cbotulinum" in patho_file or "Cbotulinum" in decoy_file:
            sim_key = "Clostridium"
        elif "Ecoli" in patho_file or "Ecoli" in decoy_file:
            sim_key = "Escherichia"
        elif "Ypestis" in patho_file or "Ypestis" in decoy_file:
            sim_key = "Yersinia"
        else: # F.tularensis and all others default to Francisella-free
            sim_key = "Francisella"

        sim_file_R1, sim_file_R2 = sim_files[sim_key]
        
        # Extract details for the output file names
        output_prefix = "WW_" + patho_file.split('_', 1)[1].rsplit('_R1', 1)[0]
        
        # --- R1 Processing ---
        print(f"\n--- Starting R1 Process for {output_prefix} ---")
        
        # 1. Process and reorder Pathogen R1 sequences (writes to uncompressed temp file)
        patho_temp_R1 = os.path.join(output_dir, f"temp_{output_prefix}_patho_R1.fastq")
        try:
            processed_patho_R1 = process_and_reorder_fastq(patho_path_R1, numbers_to_keep, patho_temp_R1)
        except Exception:
            print(f"Skipping pair {output_prefix} due to R1 processing error.")
            continue # Move to the next pair

        # 2. Process and reorder Decoy R1 sequences (writes to uncompressed temp file)
        decoy_temp_R1 = os.path.join(output_dir, f"temp_{output_prefix}_decoy_R1.fastq")
        try:
            processed_decoy_R1 = process_and_reorder_fastq(decoy_path_R1, numbers_to_keep, decoy_temp_R1)
        except Exception:
            print(f"Skipping pair {output_prefix} due to R1 processing error.")
            # Clean up the successfully processed file
            os.remove(patho_temp_R1)
            continue # Move to the next pair
            
        # The background metagenome (sim_file_R1) is used AS IS and is gzipped, 
        # so we need to decompress it to a temp file for interleaving.
        sim_temp_R1 = os.path.join(output_dir, f"temp_{output_prefix}_sim_R1.fastq")
        print(f"Decompressing {os.path.basename(sim_file_R1)} to {os.path.basename(sim_temp_R1)}")
        with gzip.open(sim_file_R1, 'rt') as f_in, open(sim_temp_R1, 'w') as f_out:
            f_out.write(f_in.read())


        # 3. Interleave R1 files (temp files are now uncompressed)
        output_file_R1 = os.path.join(output_dir, f"{output_prefix}_R1.fastq.gz")
        interleave_fastqs(sim_temp_R1, processed_patho_R1, processed_decoy_R1, output_file_R1)
        
        # --- R2 Processing ---
        print(f"\n--- Starting R2 Process for {output_prefix} ---")
        
        # The R2 files for Pathogen and Decoy correspond directly to R1 reads,
        # so we use the R1 logic's 'numbers_to_keep' which is based on the R1 files.
        patho_path_R2 = patho_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')
        decoy_path_R2 = decoy_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')
        
        # 1. Process and reorder Pathogen R2 sequences
        patho_temp_R2 = os.path.join(output_dir, f"temp_{output_prefix}_patho_R2.fastq")
        try:
            processed_patho_R2 = process_and_reorder_fastq(patho_path_R2, numbers_to_keep, patho_temp_R2)
        except Exception:
            print(f"Skipping R2 for pair {output_prefix} due to R2 processing error.")
            # Clean up R1 temp files
            os.remove(patho_temp_R1); os.remove(decoy_temp_R1); os.remove(sim_temp_R1)
            continue
            
        # 2. Process and reorder Decoy R2 sequences
        decoy_temp_R2 = os.path.join(output_dir, f"temp_{output_prefix}_decoy_R2.fastq")
        try:
            processed_decoy_R2 = process_and_reorder_fastq(decoy_path_R2, numbers_to_keep, decoy_temp_R2)
        except Exception:
            print(f"Skipping R2 for pair {output_prefix} due to R2 processing error.")
            # Clean up R1 temp files
            os.remove(patho_temp_R1); os.remove(decoy_temp_R1); os.remove(sim_temp_R1)
            # Clean up R2 temp files
            os.remove(patho_temp_R2)
            continue

        # Decompress Simulated R2 file
        sim_temp_R2 = os.path.join(output_dir, f"temp_{output_prefix}_sim_R2.fastq")
        print(f"Decompressing {os.path.basename(sim_file_R2)} to {os.path.basename(sim_temp_R2)}")
        with gzip.open(sim_file_R2, 'rt') as f_in, open(sim_temp_R2, 'w') as f_out:
            f_out.write(f_in.read())

        # 3. Interleave R2 files
        output_file_R2 = os.path.join(output_dir, f"{output_prefix}_R2.fastq.gz")
        interleave_fastqs(sim_temp_R2, processed_patho_R2, processed_decoy_R2, output_file_R2)
        
        # --- Cleanup ---
        print(f"\n--- Cleaning up temporary files for {output_prefix} ---")
        os.remove(patho_temp_R1)
        os.remove(decoy_temp_R1)
        os.remove(sim_temp_R1)
        os.remove(patho_temp_R2)
        os.remove(decoy_temp_R2)
        os.remove(sim_temp_R2)
        print("Cleanup complete.")


if __name__ == '__main__':
    main()
    print("\n\nInterleaving and Reordering complete. Final output files are in:", output_dir)

####################################################################################################################################################################

    #!/usr/bin/env python3

    import os
    from Bio import SeqIO
    import itertools
    import gzip
    import sys

    simulated_dir = "${params.bg_dir}"
    patho_decoy_dir = "${params.patho_decoy_dir}"
    output_dir = "${params.output_dir}"

    # create the function that interleaves the background metagenome with the pathogen and decoy files
    def interleave_fastqs(file1, file2, file3, output_file):
        # Open the output file with gzip.open so that the output is in gz compressed format
        with gzip.open(output_file, 'wt') as out_handle:
            try:
                records1 = SeqIO.parse(gzip.open(file1, 'rt'), 'fastq')     # Should the fastq be fastq.gz?
                records2 = SeqIO.parse(gzip.open(file2, 'rt'), 'fastq')
                records3 = SeqIO.parse(gzip.open(file3, 'rt'), 'fastq')
                # Make sure all files are processed regardless of length
                for r1, r2, r3 in itertools.zip_longest(records1, records2, records3):
                    if r1:
                        SeqIO.write(r1, out_handle, 'fastq')
                    if r2:
                        SeqIO.write(r2, out_handle, 'fastq')
                    if r3:
                        SeqIO.write(r3, out_handle, 'fastq')
            # Error messages
            except FileNotFoundError as e:
                print(f"Error: One of the input files was not found: {e}", file=sys.stderr)
                sys.exit(1) # Exit with an error code
            except Exception as e:
                print(f"An unexpected error occurred during file processing: {e}", file=sys.stderr)
                sys.exit(1) # Exit with an error code


    # Check that output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get the pathogen and decoy files
    all_files = os.listdir(patho_decoy_dir)
    patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq.gz")])
    decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq.gz")])

    # Get the simulated wastewater files
    sim_woBacillus_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoBacillus_R1.fastq.gz")
    sim_woBacillus_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoBacillus_R2.fastq.gz")
    sim_woClostridium_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoClostridium_R1.fastq.gz")
    sim_woClostridium_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoClostridium_R2.fastq.gz")
    sim_woEscherichia_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoEscherichia_R1.fastq.gz")
    sim_woEscherichia_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoEscherichia_R2.fastq.gz")
    sim_woFrancisella_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoFrancisella_R1.fastq.gz")
    sim_woFrancisella_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoFrancisella_R2.fastq.gz")
    sim_woYersinia_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoYersinia_R1.fastq.gz")
    sim_woYersinia_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoYersinia_R2.fastq.gz")

    # Process each pair of pathogen and decoy files
    for patho_file, decoy_file in zip(patho_files, decoy_files):
        patho_path = os.path.join(patho_decoy_dir, patho_file)
        decoy_path = os.path.join(patho_decoy_dir, decoy_file)
        
        # Extract details in the file names for the output file names
        output_prefix = "WW_" + patho_file.split('_', 1)[1].rsplit('_R1', 1)[0]
        
        # Select the correct background file for the pathogen/decoy pair
        # For the B.anthracis pathogen/decoy files
        if "Banthracis" in patho_file or "Banthracis" in decoy_file:
            sim_file_R1 = sim_woBacillus_R1
            sim_file_R2 = sim_woBacillus_R2
        # For the C.botulinum pathogen/decoy pair
        elif "Cbotulinum" in patho_file or "Cbotulinum" in decoy_file:
            sim_file_R1 = sim_woClostridium_R1
            sim_file_R2 = sim_woClostridium_R2
        # For the E.coli pathogen/decoy pair
        elif "Ecoli" in patho_file or "Ecoli" in decoy_file:
            sim_file_R1 = sim_woEscherichia_R1
            sim_file_R2 = sim_woEscherichia_R2
        # For the Y.pestis pathogen/decoy pair
        elif "Ypestis" in patho_file or "Ypestis" in decoy_file:
            sim_file_R1 = sim_woYersinia_R1
            sim_file_R2 = sim_woYersinia_R2
        # For all other files including F.tularensis
        else:
            sim_file_R1 = sim_woFrancisella_R1
            sim_file_R2 = sim_woFrancisella_R2

        # Interleave R1 files
        output_file_R1 = os.path.join(output_dir, f"{output_prefix}_R1.fastq.gz")
        interleave_fastqs(sim_file_R1, patho_path, decoy_path, output_file_R1)
        
        # Change suffix from R1 to R2
        patho_file_R2 = patho_path.replace('_R1.fastq.gz', '_R2.fastq.gz')
        decoy_file_R2 = decoy_path.replace('_R1.fastq.gz', '_R2.fastq.gz')
        # Interleave R2 files
        output_file_R2 = os.path.join(output_dir, f"{output_prefix}_R2.fastq.gz")
        interleave_fastqs(sim_file_R2, patho_file_R2, decoy_file_R2, output_file_R2)

    print("Interleaving complete. Output files are in:", output_dir)

######################################################################################################################################################################

#!/usr/bin/env python3

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
output_dir = "${params.output_dir}"

# create the function that interleaves the background metagenome with the pathogen and decoy files into a single file for each pathogen/decoy pair
def interleave_fastqs(file1, file2, file3, output_path):
    with gzip.open(output_path, 'wt') as out_handle:
        try:
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
        except FileNotFoundError as e:
            print(f"Error: Input file not found: {e}", file=sys.stderr)
            sys.exit(1) # Exit with an error code
        except Exception as e:
            print(f"Error during processing: {e}", file=sys.stderr)
            sys.exit(1) # Exit with an error code

# make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# create a look up table that maps pathogens to their background "without" files
bg_map = {
    "Banthracis": "Simulated_WastewaterwoBacillus",
    "Cbotulinum": "Simulated_WastewaterwoClostridium",
    "Ecoli":      "Simulated_WastewaterwoEscherichia",
    "Ypestis":    "Simulated_WastewaterwoYersinia",
    "default":    "Simulated_WastewaterwoFrancisella"
}

# sort and pair up the pathogen and decoy pairs
all_files = os.listdir(patho_decoy_dir)
patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq.gz")])
decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq.gz")])

# process each pathogen decoy pair
for patho_file, decoy_file in zip(patho_files, decoy_files):
    patho_path_R1 = os.path.join(patho_decoy_dir, patho_file)
    decoy_path_R1 = os.path.join(patho_decoy_dir, decoy_file)
    
    # search loop to determine which background file to use
    bg_prefix = bg_map["default"]
    for key in bg_map:
        if key in patho_file:
            bg_prefix = bg_map[key]
            break

    # construct background and R2 paths
    sim_file_R1 = os.path.join(simulated_dir, f"{bg_prefix}_R1.fastq.gz")
    sim_file_R2 = os.path.join(simulated_dir, f"{bg_prefix}_R2.fastq.gz")
    patho_path_R2 = patho_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')
    decoy_path_R2 = decoy_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')

    output_prefix = "WW_" + patho_file.split('_', 1)[1].rsplit('_R1', 1)[0]

    # --- Use Temporary Files ---
    # NamedTemporaryFile allows other processes to find the file by name
    with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as tmp_R1, \
         tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as tmp_R2:
        
        print(f"Processing {output_prefix} via temps: {tmp_R1.name}")

        # Interleave into the temporary files
        interleave_fastqs(sim_file_R1, patho_path_R1, decoy_path_R1, tmp_R1.name)
        interleave_fastqs(sim_file_R2, patho_path_R2, decoy_path_R2, tmp_R2.name)

        # Move the temporary files to the final destination
        final_R1 = os.path.join(output_dir, f"{output_prefix}_R1.fastq.gz")
        final_R2 = os.path.join(output_dir, f"{output_prefix}_R2.fastq.gz")
        
        shutil.move(tmp_R1.name, final_R1)
        shutil.move(tmp_R2.name, final_R2)

print(f"Interleaving complete. Files moved to: {output_dir}")


#######################################################################################################################################################################################

#!/usr/bin/env python3

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
output_dir = "${params.output_dir}"

# create the function that interleaves the background metagenome with the pathogen and decoy files into a single file for each pathogen/decoy pair
def interleave_fastqs(file1, file2, file3, output_path):
    with gzip.open(output_path, 'wt') as out_handle:
        try:
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
        except FileNotFoundError as e:
            print(f"Error: Input file not found: {e}", file=sys.stderr)
            sys.exit(1) # Exit with an error code
        except Exception as e:
            print(f"Error during processing: {e}", file=sys.stderr)
            sys.exit(1) # Exit with an error code

# make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# create a look up table that maps pathogens to their background "without" files
bg_map = {
    "Banthracis": "Simulated_WastewaterwoBacillus",
    "Cbotulinum": "Simulated_WastewaterwoClostridium",
    "Ecoli":      "Simulated_WastewaterwoEscherichia",
    "Ypestis":    "Simulated_WastewaterwoYersinia",
    "default":    "Simulated_WastewaterwoFrancisella"
}

# sort and pair up the pathogen and decoy pairs
all_files = os.listdir(patho_decoy_dir)
patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq.gz")])
decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq.gz")])

    # store the paths to the temporary files to pass to the next function
    interleaved_results = []

    for patho_file, decoy_file in zip(patho_files, decoy_files):
        # search loop to determine which background file to use
        bg_prefix = bg_map["default"]
        for key in bg_map:
            if key in patho_file:
                bg_prefix = bg_map[key]
                break

        # setup input paths (background, pathogen an decoy)
        sim_file_R1 = os.path.join(simulated_dir, f"{bg_prefix}_R1.fastq.gz")
        sim_file_R2 = os.path.join(simulated_dir, f"{bg_prefix}_R2.fastq.gz")
        patho_path_R1 = os.path.join(patho_decoy_dir, patho_file)
        patho_path_R2 = patho_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')
        decoy_path_R1 = os.path.join(patho_decoy_dir, decoy_file)
        decoy_path_R2 = decoy_path_R1.replace('_R1.fastq.gz', '_R2.fastq.gz')

        # --- Create Identifiable Temporary Files ---
        # prefix: makes it recognizable (WW_Banthracis_randomID)
        # suffix: helps other tools identify the file type
        # delete=False: keeps it for the next function
        tmp_R1 = tempfile.NamedTemporaryFile(prefix=f"{output_prefix}_", suffix="_R1.fastq.gz", delete=False)
        tmp_R2 = tempfile.NamedTemporaryFile(prefix=f"{output_prefix}_", suffix="_R2.fastq.gz", delete=False)
        
        tmp_R1.close()
        tmp_R2.close()

        print(f"Creating temp file: {tmp_R1.name}")

        interleave_fastqs(sim_file_R1, patho_path_R1, decoy_path_R1, tmp_R1.name)
        interleave_fastqs(sim_file_R2, patho_path_R2, decoy_path_R2, tmp_R2.name)

        interleaved_results.append({
            "sample_id": output_prefix,
            "r1_path": tmp_R1.name,
            "r2_path": tmp_R2.name
        })

    return interleaved_results

def next_analysis_step(data_list):
    """Example function receiving the temp files."""
    for entry in data_list:
        print(f"Next step processing: {entry['sample_id']} at {entry['r1_path']}")

if __name__ == "__main__":
    results = process_all_files()
    next_analysis_step(results)

##########################################################################################################################################################################

import os
import gzip
import tempfile
from Bio import SeqIO

# Function to read in the numbers for the sequences
def read_number_list(number_list_file):
    """
    read the input number list from the given text file and return a list of integers
    """
    try:
        with open(number_list_file, 'r') as f:
            # read all lines, strip whitespace, and convert to integer
            # handle files where numbers might be on separate lines or space/comma separated
            numbers = [int(num) for line in f for num in line.strip().split() if num.isdigit()]
        return numbers
    except FileNotFoundError:
        print(f"Error: {number_list_file} does not exist", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: cannot read {e}", file=sys.stderr)
        sys.exit(1)

# --- Function to shuffle reads  based on ---
def shuffle_and_extract_reads(data_list, number_list):
    shuffled_results = []

    for entry in data_list:
        print(f"Shuffling reads for {entry['sample_id']}...")

        # Create new temp files for the shuffled output
        shuffled_R1 = tempfile.NamedTemporaryFile(prefix=f"shuffled_{entry['sample_id']}_", suffix="_R1.fastq.gz", delete=False)
        shuffled_R2 = tempfile.NamedTemporaryFile(prefix=f"shuffled_{entry['sample_id']}_", suffix="_R2.fastq.gz", delete=False)
        shuffled_R1.close()
        shuffled_R2.close()

        # load the interleaved reads into a temporary index or dictionary
        read_map_R1 = {}
        read_map_R2 = {}

        with gzip.open(entry['r1_path'], 'rt') as f1, \
             gzip.open(entry['r2_path'], 'rt') as f2:
            
            # enumerate starts at 1 to match 'numbers' list logic
            for i, (rec1, rec2) in enumerate(zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')), 1):
                read_map_R1[i] = rec1
                read_map_R2[i] = rec2

        # write out the reads in the order specified by number_list
        with gzip.open(shuffled_R1.name, 'wt') as out1, \
             gzip.open(shuffled_R2.name, 'wt') as out2:
            
            for num in number_list:
                # skip if the number is higher than the available reads
                if num in read_map_R1:
                    SeqIO.write(read_map_R1[num], out1, 'fastq')
                    SeqIO.write(read_map_R2[num], out2, 'fastq')
                else:
                    # skip indices that exceed the file length as requested
                    continue

        shuffled_results.append({
            "sample_id": f"shuffled_{entry['sample_id']}",
            "r1_path": shuffled_R1.name,
            "r2_path": shuffled_R2.name
        })
        
        # clean up the original interleaved temp files to save space
        os.remove(entry['r1_path'])
        os.remove(entry['r2_path'])

    return shuffled_results

# --- How to use it in your main block ---
if __name__ == "__main__":
    # 1. Generate the interleaved files
    initial_results = process_all_files()
    
    # 2. Read your list of desired read positions
    number_list_file = "${params.numbers}" # Path from your config
    numbers = read_number_list(number_list_file)
    
    # 3. Shuffle them!
    final_shuffled_data = shuffle_and_extract_reads(initial_results, numbers)
    
    # Now final_shuffled_data contains the paths to your ordered reads
    print("Shuffle complete.")

############################################################################################################################################################################

#!/usr/bin/env python3

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
        try:
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
        except FileNotFoundError as e:
            print(f"Error: Input file not found: {e}", file=sys.stderr)
            sys.exit(1) # Exit with an error code
        except Exception as e:
            print(f"Error during processing: {e}", file=sys.stderr)
            sys.exit(1) # Exit with an error code

# --- Function to read in the numbers for the sequences ---
def read_number_list(number_list_file):
    try:
        with open(number_list_file, 'r') as f:
            # read all lines, strip whitespace, and convert to integer
            # handle files where numbers might be on separate lines or space/comma separated
            numbers = [int(num) for line in f for num in line.strip().split() if num.isdigit()]
        return numbers
    except FileNotFoundError:
        print(f"Error: {number_list_file} does not exist", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: cannot read {e}", file=sys.stderr)
        sys.exit(1)

# --- Function to determine which files to interleave together ---
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

# --- Function to shuffle the interleaved files ---
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

# --- Export shuffled files ---
def export_final_files(data_list, destination):
    # moves files from /tmp to permanent storage with final file names
    os.makedirs(destination, exist_ok=True)
    for entry in data_list:
        final_R1 = os.path.join(destination, f"spiked_{entry['sample_id']}_R1.fastq.gz")
        final_R2 = os.path.join(destination, f"spiked_{entry['sample_id']}_R2.fastq.gz")
        
        shutil.move(entry['r1_path'], final_R1)
        shutil.move(entry['r2_path'], final_R2)
        print(f"Exported: {final_R1}")

# ---ALTERNATIVELY (I PREFER) ---
# --- Export Shuffled Files with Subdirectory Logic ---
def export_final_files(data_list, destination):
    # moves files from /tmp to permanent storage with folder-specific names
    for entry in data_list:
        sample_id = entry['sample_id']
        
        # Determine the correct subdirectory based on the name
        sub_folder = ""
        if "_f1_" in sample_id.lower():
            sub_folder = "F1"
        elif "_f2_" in sample_id.lower():
            sub_folder = "F2"
        elif "_f3_" in sample_id.lower():
            sub_folder = "F3"
        
        # Construct the full target path
        target_dir = os.path.join(destination, sub_folder)
        
        # Create the subdirectory (and parent destination) if they don't exist
        os.makedirs(target_dir, exist_ok=True)
        
        final_R1 = os.path.join(target_dir, f"spiked_{sample_id}_R1.fastq.gz")
        final_R2 = os.path.join(target_dir, f"spiked_{sample_id}_R2.fastq.gz")
        
        # Move the files from temp to the organized subfolders
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
    
    print("\nPipeline finished successfully.")
