#!/usr/bin/env python3

import subprocess
import os
import glob
import sys

# Input directory
input_dir = "<YOUR_INPUT_DIRECTORY>"

# Output directory
output_dir = os.path.join(input_dir, "COUNTS_OUT")
os.makedirs(output_dir, exist_ok=True)

# GTF file
gtf_file = "<YOUR.PATH.TO.GTF.GENOME>"

# Get list of BAM files
file_list = sorted(glob.glob(os.path.join(input_dir, "*Aligned.sortedByCoord.out.bam")))

# SLURM task ID
try:
    array_task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID"))
    bam_file = file_list[array_task_id]
except (TypeError, IndexError):
    print("Error: Invalid or missing SLURM_ARRAY_TASK_ID.")
    sys.exit(1)

# Output file
output_file = os.path.splitext(os.path.basename(bam_file))[0] + "_counts.txt"
output_path = os.path.join(output_dir, output_file)

# HTSeq command (GTF uses exons + gene_id fields)
htseq_command = [
    "htseq-count",
    "-f", "bam",
    "-r", "pos",
    "-s", "no",              
    "-t", "exon",            
    "-i", "gene_id",         
    bam_file,
    gtf_file
]

# Run HTSeq
try:
    with open(output_path, "w") as output_stream:
        result = subprocess.run(htseq_command, stdout=output_stream, check=True, stderr=subprocess.PIPE)
    print(f"HTSeq counting completed successfully for {bam_file}")
except subprocess.CalledProcessError as e:
    print(f"Error running HTSeq on {bam_file}: {e}")
    print(f"Stderr output:\n{e.stderr.decode()}")
