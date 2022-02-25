

#job_directory=$PWD/.job

job_directory=/projects/b1059/projects/Ryan/protein_structure/af_binary_pred/bicoid/

files=$job_directory*.fa
#echo $files

for filename in $files; do 
    #echo $filename
    file=$(basename $filename) #Get just the file
    #echo $file
    name=${file%.*} #Remove extension
    #echo $name

    input_fasta=$filename

    #echo $job_file

    #Create a submission file
    cpu_file="${job_directory}${name}_cpu.sh"
    echo "#!/bin/bash
#SBATCH --account=b1042  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=genomics  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=12 ## how many cpus or processors do you need on each computer
#SBATCH --time=3:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem=80G ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=$name  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --output=$name.log ## standard out and standard error goes to this file

#########################################################################
### PLEASE NOTE:                                                      ###
### The above CPU, Memory, and GPU resources have been selected based ###
### on the computing resources that alphafold was tested on           ###
### which can be found here:                                          ###
### https://github.com/deepmind/alphafold#running-alphafold)          ###
### It is likely that you do not have to change anything above        ###
### besides your allocation, and email (if you want to be emailed).   ###
#########################################################################

module purge
module load alphafold/2.1.1-only-msas-flag-addition

# template
# real example
alphafold-multimer --fasta_paths=$input_fasta\
    --output_dir=$job_directory \
    --max_template_date=2021-11-01 \
    --model_preset=multimer \
    --db_preset=full_dbs  \
    --is_prokaryote_list=false \
    --max_template_date=2021-11-01 \
    --only_msas=true " > $cpu_file


gpu_file="${job_directory}${name}_gpu.sh"
echo "#!/bin/bash
#SBATCH --account=b1042  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=genomics-gpu  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on each computer
#SBATCH --gres=gpu:a100:1  ## type of GPU requested, and number of GPU cards to run on
#SBATCH --time=1:30:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem=85G ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=$name  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --output=$name.log ## standard out and standard error goes to this file

#########################################################################
### PLEASE NOTE:                                                      ###
### The above CPU, Memory, and GPU resources have been selected based ###
### on the computing resources that alphafold was tested on           ###
### which can be found here:                                          ###
### https://github.com/deepmind/alphafold#running-alphafold)          ###
### It is likely that you do not have to change anything above        ###
### besides your allocation, and email (if you want to be emailed).   ###
#########################################################################

module purge
module load alphafold/2.1.1-only-msas-flag-addition

# template
# real example
alphafold-multimer --fasta_paths=$input_fasta \
    --max_template_date=2021-11-01 \
    --model_preset=multimer \
    --db_preset=full_dbs  \
    --use_precomputed_msas=true \
    --output_dir $job_directory \
    --is_prokaryote_list=false " > $gpu_file

#Submite the actual jobs
cpu_job=($(sbatch $cpu_file))

echo "cpu_job ${cpu_job[-1]}" >> slurm_ids

gpu_job=($(sbatch --dependency=afterok:${cpu_job[-1]} $gpu_file))

echo "gpu_job ${gpu_job[-1]}" >> slurm_ids

done



