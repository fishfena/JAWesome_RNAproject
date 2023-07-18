#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

def run_fastqc(folder_path):
    # Get a list of all .fastq.gz files in the folder
    files = [f for f in os.listdir(folder_path) if f.endswith('.fastq.gz')]

    # Loop through each file and run FastQC
    for file in files:
        # Build the full path to the input file
        input_file = os.path.join(folder_path, file)

 # Trim adapters using TrimGalore
        trim_output_prefix = os.path.splitext(file)[0]
        trim_output_folder = os.path.join(folder_path, trim_output_prefix + "_trimmed")
        trim_output_file = os.path.join(trim_output_folder, trim_output_prefix + "_trimmed.fq.gz")

        trim_cmd = ['trim_galore', '--output_dir', trim_output_folder, '--gzip', input_file]

        # Run TrimGalore using subprocess
        subprocess.run(trim_cmd)

        print(f"TrimGalore completed for {file}")

# Example usage
folder_path = '/path/to/your/folder/in/the/cluster/or/your/local/machine/'
run_fastqc(folder_path)
        
