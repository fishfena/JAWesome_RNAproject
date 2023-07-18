#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 13:50:22 2023

@author: mfp
"""

import os

# Set the paths
genome_dir = "/genome/dir/path/"
trimmed_seq_dir = "/trimmed/sequences/"

# Set the number of threads for STAR
threads = 4

# Loop through the fastq files in the trimmed_seq_dir
for fastq_file in os.listdir(trimmed_seq_dir):
    if fastq_file.endswith("_trimmed.fq.gz") and not fastq_file.endswith(".Z"):
        # Extract the sample name from the file name
        sample_name = fastq_file.replace("_R1_001_trimmed.fq.gz", "")
        
        # Run STAR alignment
        star_command = (
            f"STAR --runThreadN {threads} "
            f"--genomeDir {genome_dir} "
            f"--readFilesCommand gzcat "
            f"--readFilesIn {os.path.join(trimmed_seq_dir, fastq_file)} "
            f"--outFileNamePrefix {os.path.join(trimmed_seq_dir, sample_name + '_')} "
            "--outSAMtype SAM"
        )
        os.system(star_command)
        
        # Convert SAM to BAM
        sam_to_bam_command = (
            f"samtools view -Sb {os.path.join(trimmed_seq_dir, sample_name + '_Aligned.out.sam')} "
            f"> {os.path.join(trimmed_seq_dir, sample_name + '.bam')}"
        )
        os.system(sam_to_bam_command)
        
        # Sort BAM
        sort_bam_command = (
            f"samtools sort {os.path.join(trimmed_seq_dir, sample_name + '.bam')} "
            f"-o {os.path.join(trimmed_seq_dir, sample_name + '.sort.bam')}"
        )
        os.system(sort_bam_command)
        
        # Index BAM
        index_bam_command = f"samtools index {os.path.join(trimmed_seq_dir, sample_name + '.sort.bam')}"
        os.system(index_bam_command)
