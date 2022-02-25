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
### Inputs 

