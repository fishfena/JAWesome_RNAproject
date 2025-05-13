#!/bin/bash
#SBATCH --job-name=<your.python.file.name>
#SBATCH --account=<your.account>
#SBATCH --partition=<your.partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=<your.python.file.name>_%A_%a.out
#SBATCH --error=<your.pythong.file.name>%A_%a.err
#SBATCH --array=0-total bam files #(i.e., 0-19 if 20 bam files) 

module load anaconda3/2024.02-1
module load bio/samtools/1.17-gcc-11.4.0
source activate htseq_env

python <your.python.file.name>.py ${SLURM_ARRAY_TASK_ID}
