#!/usr/bin/env nextflow

params.input_csv = "/experiment/chain_key.csv"
params.fasta = "some/path/monomer.fa"



process processInputs {
    
    
    """
    python3 scripts/process_inputs.py -p \$parsing_file -fm \$mon_fasta \$-o
    
    """

}

process submitJobs {

    input:
    path 'jobdir' from 

    """
    bash parallel_submission.sh
    
    """

}