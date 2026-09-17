#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import Helper
// import CheckParams (LMAS lib for tutorial)


/*
 * Help message
*/
// Help message (see print_help section in LMAS_main/lib/Helper.groovy for help lol)
params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

/*
 * Checksum
*/
process verify_files {
    tag "${all_files.name}"
    
    input:
    path all_files
    path checksum_sha

    script:
    """
    #!/usr/bin/env python3
    import hashlib
    import csv
    import sys

    filename = "${all_files.name}"
    expected_hash = ""

    #  look up the hash in the CSV
    with open("${checksum_sha}", mode='r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == filename:
                expected_hash = row[1].strip().lower()
                break

    if not expected_hash:
        print(f"ERROR: {filename} not found in CSV.")
        sys.exit(1)

    #  calculate the SHA256 hash of the file
    sha256_hash = hashlib.sha256()
    with open("${all_files}", "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    
    calculated_hash = sha256_hash.hexdigest().lower()

    # compare the calculated hash to the list hash
    if calculated_hash == expected_hash:
        print(f"SUCCESS: {filename} matches.")
    else:
        print(f"ERROR: Mismatch for {filename}!")
        print(f"Expected: {expected_hash}")
        print(f"Found:    {calculated_hash}")
        sys.exit(1)
    """
}

/*
 * Merge and shuffle process
*/
process shuffle {
    label 'SCS'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path(ww_files)
    path(pd_files)
    val(rand_num_list)

    output:
    path('spiked_*.fastq.gz')

    script:
    """
    #!/usr/bin/env python3
    """
}

/*
 * Run the main workflow
*/
workflow{
    // Main parameters
    wastewater_pairs = Channel.fromFilePairs("${params.bg_dir}/Simulated_Wastewater*_{R1,R2}.fastq.gz").ifEmpty {exit 1, "No such fastq files in path:'${params.bg_dir}'/Simulated_Wastewater..."}
    patho_decoy_pairs = Channel.fromFilePairs("${params.patho_decoy_dir}/**/{patho,decoy}_*{R1,R2}.fastq.gz").ifEmpty {exit 1, "No such fastq files in path:'${params.patho_decoy_dir}'patho... or decoy..."}
    number_list_file = Channel.fromPath(params.numbers).ifEmpty {exit 1, "No file found with pattern:'${params.numbers}'"}
    checker =  Channel.fromPath(params.check_sum_list).splitCsv(header:true).ifEmpty {exit 1, "No file found with pattern:'${params.check_sum_list}'"}

    // Checksum step
    verify_files()

    // Merge and shuffle steps
    shuffle(wastewater_pairs, patho_decoy_pairs, number_list_file)
}
workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}
workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
