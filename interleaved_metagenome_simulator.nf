#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// this prevents a warning of undefined parameter
//params.help             = false

/*
 * Generate simulated reads for the microbes in the metagenome
 * Channel in the csv file containing genomes, read counts and location of the microbes used for the metagenome
*/
Channel
  .fromPath(params.bkgenv_infile, checkIfExists:true)
  .splitCsv(header:true)
  .set {metagenome_ch}
metagenome_ch.view()

/*
 * Channel in the csv file containing genomes, read counts and location of the microbes used for the pathogen and decoy bacteria
*/
Channel
  .fromPath(params.spikein_file, checkIfExists:true)
  .splitCsv(header:true)
  .set {spikein_ch}
spikein_ch.view()

/*
 * Define the number of iterations for the pathogen and non-pathogen
*/ 
iterate = Channel.value(1..params.iteration)

/*
 * Process for generating the reads for the metagenome
*/
process simulate_environment {
    storeDir ${params.interleaved_sim_meta}
    input:
    tuple val(Name), val(refseq_id), val(Genome_Type), val(single_read_count), val(Read_Count), val(Filename)
    output:
    path("mg_${Microbial_Name}_*.fastq"), emit: microbial_fastqs
    script:
    """
    iss generate ${Genome_Type} ${Filename} --n_reads ${Read_Count} --model miseq \
    --abundance uniform --mode kde --debug --seed 125678 --output mg_${Microbial_Name}
    """
}

/*
 * Merge the reads for the environmental microbes to create the simulated metagenome

process background_metagenome {
//    storeDir ${params.interleaved_sim_meta}
    input:
    path microbial_fastqs
    output:
    file"${params.mgenv_outfile_R1}", emit: forward_metagenome
    file"${params.mgenv_outfile_R2}", emit: reverse_metagenome
    script:
    """
    cat *_R1.fastq > ${params.mgenv_outfile_R1}
    cat *_R2.fastq > ${params.mgenv_outfile_R2}
    """
}
*/

/*
 * Process for generating the reads for the pathogen and non-pathogen
*/
process simulate_pathogen {
  input:
  tuple val(Patho_Name), path(Patho_Filename), val(Patho_Read_Count)
  each(iterate)
  file (patho_abundanceFile)
  file (decoy_abundanceFile)
  output:
  path("${Patho_Name}_${iter}*.fastq"), emit: pathogen_reads

  script:
  def file = (Patho_Name)
  def abundanceFile = ''
  
  if (file.startsWith("patho")) {
      abundanceFile = (patho_abundanceFile)
  } else if (file.startsWith("decoy")) {
      abundanceFile = (decoy_abundanceFile)
      println "Invalid input file name"
  }
  
  println "abundance file is: $abundanceFile"

    """
    iss generate --genomes ${Patho_Filename} --n_reads ${Patho_Read_Count} --model miseq \
    --abundance_file ${abundanceFile} --mode kde --debug --seed ${iterate} --output ${Patho_Name}_${iterate}
    """
}

/*
 * Option A: Process to interleave the reads of the metagenome
 * More favoured strategy for full generation and distribution interleaving.
 * Contains config file contents and .nf file content
*/

params{
    input_dir = "path/to/bgmetagenome/fastqs"
    patho_decoy_dir = "path/to/patho/and/decoy/fastqs"
    output_dir = "path/to/output/folder"
    num_files = 10 // Default number of files to process, can be overridden
}

// Define input channels
Channel
    .fromFilePairs("${params.input_dir}/Simulated_Wastewater*_{R1,R2}.fastq")
    .set { wastewater_pairs }

Channel
    .fromFilePairs("${params.patho_decoy_dir}/{patho,decoy}_*_{R1,R2}.fastq")
    .take( params.num_files )
    .set { patho_decoy_pairs }

process interleaveFiles {
    storeDir "${params.output_dir}"

    input:
    tuple val(ww_name), path(ww_files) from wastewater_pairs
    tuple val(pd_name), path(pd_files) from patho_decoy_pairs

    output:
    path "WW_*.fastq"

    script:
    """
    #!/usr/bin/env python3

    import os
    from Bio import SeqIO
    import itertools
    import glob

    input_dir = "/home/tj/Downloads/Benchmarking/Generator/test/Dummy_fastq"
    patho_decoy_dir = "/home/tj/Downloads/Benchmarking/Generator/test/Dummy_pathodecoy"
    output_dir = "/home/tj/Downloads/Benchmarking/Generator/test/Dummy_out"

    def interleave_fastqs(file1, file2, file3, output_file):
        with open(output_file, 'w') as out_handle:
            try:
                records1 = SeqIO.parse(file1, 'fastq')
                records2 = SeqIO.parse(file2, 'fastq')
                records3 = SeqIO.parse(file3, 'fastq')
                for r1, r2, r3 in itertools.zip_longest(records1, records2, records3):
                    if r1:
                        SeqIO.write(r1, out_handle, 'fastq')
                    if r2:
                        SeqIO.write(r2, out_handle, 'fastq')
                    if r3:
                        SeqIO.write(r3, out_handle, 'fastq')
            except Exception as e:
                print(f"Error processing files: {e}")

    # Paths to the directories
    patho_decoy_dir = patho_decoy_dir
    simulated_dir = input_dir
    output_dir = output_dir

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get the patho and decoy files
    all_files = os.listdir(patho_decoy_dir)
    patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq")])
    decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq")])

    # Get the simulated files
    sim_wEcoli_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwEcoli_R1.fastq")
    sim_wEcoli_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwEcoli_R2.fastq")
    sim_woEcoli_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoEcoli_R1.fastq")
    sim_woEcoli_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoEcoli_R2.fastq")

    # Process each pair of patho and decoy files
    for patho_file, decoy_file in zip(patho_files, decoy_files):
        patho_path = os.path.join(patho_decoy_dir, patho_file)
        decoy_path = os.path.join(patho_decoy_dir, decoy_file)
        
        output_prefix = "WW_" + patho_file.split('_', 1)[1].rsplit('_R1', 1)[0]
        
        if "Ecoli" in patho_file or "Ecoli" in decoy_file:
            sim_file_R1 = sim_wEcoli_R1
            sim_file_R2 = sim_wEcoli_R2
        else:
            sim_file_R1 = sim_woEcoli_R1
            sim_file_R2 = sim_woEcoli_R2

        # Interleave R1 files
        output_file_R1 = os.path.join(output_dir, f"{output_prefix}_R1.fastq")
        interleave_fastqs(sim_file_R1, patho_path, decoy_path, output_file_R1)
        
        # Interleave R2 files
        patho_file_R2 = patho_path.replace('_R1.fastq', '_R2.fastq')
        decoy_file_R2 = decoy_path.replace('_R1.fastq', '_R2.fastq')
        output_file_R2 = os.path.join(output_dir, f"{output_prefix}_R2.fastq")
        interleave_fastqs(sim_file_R2, patho_file_R2, decoy_file_R2, output_file_R2)

    print("Interleaving complete. Output files are in:", output_dir)
    """
}


/*
 * Option B: Process to interleave the reads of the metagenome
 * More favoured strategy for full generation and distribution interleaving.
 * Contains config file contents and .nf file content
 * for gzip fastq files
 * alternative to allow files that are empty
*/

params{
    input_dir = ""
    patho_decoy_dir = ""
    output_dir = ""
    num_files = 10 // Default number of files to process, can be overridden
}

// Create output directory
storeDir params.output_dir

// Define input channels
Channel
    .fromFilePairs("${params.patho_decoy_dir}/{patho,decoy}_*.fastq", size: 2)
    .take(params.num_files)
    .set { patho_decoy_pairs }

Channel
    .fromFilePairs("${params.simulated_dir}/Simulated_Wastewater*_R{1,2}.fastq", size: 2)
    .set { simulated_pairs }

process interleaveFiles {
    errorStrategy 'ignore'

    input:
    tuple val(pair_id), path(patho_decoy_files) from patho_decoy_pairs
    tuple val(sim_id), path(simulated_files) from simulated_pairs

    output:
    tuple val(pair_id), path("WW_*.fastq") into interleaved_output

    script:
    """
    #!/usr/bin/env python3

    import os
    from Bio import SeqIO
    import itertools
    import glob

    input_dir = "/home/tj/Downloads/Benchmarking/Generator/test/Dummy_fastq"
    patho_decoy_dir = "/home/tj/Downloads/Benchmarking/Generator/test/Dummy_pathodecoy"
    output_dir = "/home/tj/Downloads/Benchmarking/Generator/test/Dummy_out"

    def interleave_fastqs(file1, file2, file3, output_file):
        with open(output_file, 'w') as out_handle:
            try:
                records1 = SeqIO.parse(file1, 'fastq')
                records2 = SeqIO.parse(file2, 'fastq')
                records3 = SeqIO.parse(file3, 'fastq')
                for r1, r2, r3 in itertools.zip_longest(records1, records2, records3):
                    if r1:
                        SeqIO.write(r1, out_handle, 'fastq')
                    if r2:
                        SeqIO.write(r2, out_handle, 'fastq')
                    if r3:
                        SeqIO.write(r3, out_handle, 'fastq')
            except Exception as e:
                print(f"Error processing files: {e}")

    # Paths to the directories
    patho_decoy_dir = "${params.patho_decoy_dir}"
    simulated_dir = "${params.simulated_dir}"
    output_dir = "${params.output_dir}"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get the patho and decoy files
    all_files = os.listdir(patho_decoy_dir)
    patho_files = sorted([f for f in all_files if f.startswith("patho_") and f.endswith("_R1.fastq")])
    decoy_files = sorted([f for f in all_files if f.startswith("decoy_") and f.endswith("_R1.fastq")])

    # Get the simulated files
    sim_wEcoli_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwEcoli_R1.fastq")
    sim_wEcoli_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwEcoli_R2.fastq")
    sim_woEcoli_R1 = os.path.join(simulated_dir, "Simulated_WastewaterwoEcoli_R1.fastq")
    sim_woEcoli_R2 = os.path.join(simulated_dir, "Simulated_WastewaterwoEcoli_R2.fastq")

    # Process each pair of patho and decoy files
    for patho_file, decoy_file in zip(patho_files, decoy_files):
        patho_path = os.path.join(patho_decoy_dir, patho_file)
        decoy_path = os.path.join(patho_decoy_dir, decoy_file)
        
        output_prefix = "WW_" + patho_file.split('_', 1)[1].rsplit('_R1', 1)[0]
        
        if "Ecoli" in patho_file or "Ecoli" in decoy_file:
            sim_file_R1 = sim_wEcoli_R1
            sim_file_R2 = sim_wEcoli_R2
        else:
            sim_file_R1 = sim_woEcoli_R1
            sim_file_R2 = sim_woEcoli_R2

        # Interleave R1 files
        output_file_R1 = os.path.join(output_dir, f"{output_prefix}_R1.fastq")
        interleave_fastqs(sim_file_R1, patho_path, decoy_path, output_file_R1)
        
        # Interleave R2 files
        patho_file_R2 = patho_path.replace('_R1.fastq', '_R2.fastq')
        decoy_file_R2 = decoy_path.replace('_R1.fastq', '_R2.fastq')
        output_file_R2 = os.path.join(output_dir, f"{output_prefix}_R2.fastq")
        interleave_fastqs(sim_file_R2, patho_file_R2, decoy_file_R2, output_file_R2)

    print("Interleaving complete. Output files are in:", output_dir)
    """
}

interleaved_output.view()


workflow {
    interleave_fastq()
    interleaved_output.view()
}


// Run command: $nextflow run script.nf --num_files 5
// nextflow run script.nf --patho_decoy_dir /path/to/patho_decoy_files --simulated_dir /path/to/simulated_files --num_files 20


workflow {
//    take:
//    metagenome_ch

    main:
    simulate_environment(metagenome_ch)
    metagenome_interleave(simulate_environment.out.collect())

    emit:
    metagenome_fastqs = simulate_environment.out
    metagenome_interleave
}
