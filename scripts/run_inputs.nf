#!/usr/bin/env nextflow

params.pairings_csv = "/projects/b1059/projects/Ryan/protein_structure/High_Throughput_AF/data/test_output"
params.monomer_fasta = ""

pairings_csv = file(params.pairings_csv)
monomer_fasta = file(params.monomer_fasta)

process processInputs {

    input: 
    path 'out_file' from out_file
    
    output:
    stdout result into 

    """
    echo $out_file
    """
}

result.view { it.trim() } 

process makeFile {
    input

}