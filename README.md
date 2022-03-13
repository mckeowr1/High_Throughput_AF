# High_Throughput_AF

Hi - This repo contains some code to run multiple Alphafold2 Multimer on a SLURM managed HPC cluster preditions and to analyze the results! 

## MultiJob Processing

This repo tries to make programmatic use of the AlphaFold more accessible. AlphaFold is set up to infer protein structure at relatively low throughput in its current distribution. A user can predict a single structure with relative ease with basic command line knowledge and access to HPC GPU resources. Resources like collabfold make this even more accessible. There is a lot of manual information entry, and configuration setting. In both cases, before a single job can even begin. After a successful prediction, models require manual inspection of structure and confidence metrics. 

While the problems described above are of minor concern to the typical user generally interested in the structure of a single protein, these manual bottlenecks create issues for researchers that require multiple predictions. 

## Multi-Multimer Modeling

This repo (at the moment) is geared explicitly towards the prediction of multimers.

Use Cases: 
- Binary Target prediction 
- Pairwise protein interaction screening
- Multi-Configuration multi-subunit prediction

## Pipeline 

## Envs
Created with the following commands: 
Conda Version: `4.5.2`

conda create -n multiprot
conda activate multiprot
conda install -c bioconda nextflow
conda install -c anaconda pandas
conda install -c anaconda numpy
conda install -c conda-forge biopython #I did have to run this twice 

### Inputs 
To generate the files for pairwise AF multimer prediction we need a Key file to set what proteins, what region of the protein, and what pairings we will model. A template of the input can be found in the example_inputs section. 
| pair_group | Name | Monomer_id | Protein_start | Protein_end | 
| ------------- |-------------| -----| ---- | --- |
| 1 | Ku70 |   |

The pair group is a unique ID given to each entry in the file and is used to generate pairings. If you were interested in Ku70 and Ku80 modeling. You would give each a unique numeric pair ID.

We also need an input fasta with the complete sequences of all proteins that will be modeled. Ie if you are interested in just one domain of Ku70 you still need to submit the full sequence of Ku70 in the fasta file.


### Processing and AF output directory 
Processing an AF output dir requires two keyfiles 
1) Dir Key - a key file created during input processing 

| group_pairing | dir_name| key_file  |
| ------------- |-------------| -----|
| 1_2 | ku70_base.ku80_base.paxx_full | path/to/key.csv  |


2) Input Key - the key submitted during input processing

python scripts/process_output.py -p path/to/directory -d max_distance -k path/to/dir_key.csv 

