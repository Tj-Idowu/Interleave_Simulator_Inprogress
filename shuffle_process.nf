nextflow.enable.dsl=2

// Process
process shuffle{
    label 'SCS'

    input:
    path(ww_files)
    path(pd_files)
    val(rand_num_list)

    output:
    path('spiked_*.fastq.gz')

    script:
    template "shuffle_compression_script.py"
}

// Workflow
workflow shuffle_flow {
    main:
    shuffle(bg_dir,patho_decoy_dir)

    emit:
    FASTQS = shuffle.out
}