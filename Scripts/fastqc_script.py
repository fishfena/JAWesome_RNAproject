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

        # Build the full path for the output file
        output_file = os.path.join(folder_path, file + "_fastqc")

        # Construct the FastQC command
        fastqc_cmd = ['fastqc', '-o', folder_path, input_file]

        # Run FastQC using subprocess
        subprocess.run(fastqc_cmd)

        print(f"FastQC completed for {file}")

# Example usage
folder_path = '/path/to/your/folder'
run_fastqc(folder_path)



