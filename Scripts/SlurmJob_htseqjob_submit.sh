#!/bin/bash
#SBATCH --job-name=htseqCounts_execute
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=<YOUR_PARTITION>
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --output=htseqCounts_execute_%A_%a.out
#SBATCH --error=htseqCounts_execute_%A_%a.err

# Set the job array range: in this case, I have 18 (samples) or BAM files = 0–17
#SBATCH --array=0-17

# === Don't forget to Load Anaconda and activate Conda environment (you can do it on the terminal too. Check what Python v. you have and htseq. ===
module load anaconda3/2024.02-1
conda activate htseq_env

# === Run Python HTSeq-counting script ===
python htseqCounts_execute.py
