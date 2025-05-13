+		#!/usr/bin/env python3
"""	
Usage:
    python htseq_array_runner_<your.info>.py <SLURM_ARRAY_TASK_ID>
"""

import subprocess
import os
import sys

# --- Step 1: Validate input ---
if len(sys.argv) < 2:
    print("Error: SLURM_ARRAY_TASK_ID not provided.\nUsage: python htseq_array_runner_<your.info>.py <task_id>")
    sys.exit(1)

task_id = int(sys.argv[1])

# --- Step 2: Define paths ---
input_dir = "<path.to.bam.files>"
output_dir = os.path.join(input_dir, "<your.folder>")
gtf_file = "<path.to.gtf.file>Cyprinodon_variegatus.C_variegatus-1.0.110.gtf"

os.makedirs(output_dir, exist_ok=True)

# --- Step 3: Define BAM file list ---
bam_files = [
    '<file1.l02.bam>',
    '<file1.l01.bam>']

try:
    bam_filename = bam_files[task_id]
except IndexError:
    print(f"Invalid SLURM_ARRAY_TASK_ID: {task_id}")
    sys.exit(1)

bam_path = os.path.join(input_dir, bam_filename)

# --- Step 4: Check input file ---
if not os.path.exists(bam_path) or os.path.getsize(bam_path) == 0:
    print(f"ERROR: BAM file missing or empty: {bam_path}")
    sys.exit(1)

# --- Step 5: Sort BAM by name ---
name_sorted_bam = bam_path.replace(".bam", ".name_sorted.bam")
sort_command = ["samtools", "sort", "-n", "-o", name_sorted_bam, bam_path]

print(f"Sorting BAM file by name: {bam_filename}")
try:
    subprocess.run(sort_command, check=True, stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    print(f"❌ samtools sort failed:\n{e.stderr.decode()}")
    sys.exit(1)

# --- Step 6: Run HTSeq ---
new_filename = bam_filename.replace("Aligned", "neWAligned").replace(".bam", "_counts.txt")
output_path = os.path.join(output_dir, new_filename)

htseq_command = [
    "htseq-count",
    "-f", "bam",
    "-r", "name",
    "-s", "no",
    "-t", "exon",
    "-i", "gene_id",
    name_sorted_bam,
    gtf_file
]

print(f"Running HTSeq on: {name_sorted_bam}")
try:
    with open(output_path, "w") as out_f:
        subprocess.run(htseq_command, stdout=out_f, stderr=subprocess.PIPE, check=True)
    print(f"✅ Finished HTSeq → {output_path}")
except subprocess.CalledProcessError as e:
    print(f"❌ HTSeq failed:\n{e.stderr.decode()}")

# --- Optional: Cleanup ---
# os.remove(name_sorted_bam)
