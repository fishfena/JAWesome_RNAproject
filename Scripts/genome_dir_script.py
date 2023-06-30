#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#generate a new genomeParameters.txt file

import os
import subprocess

def generate_genome_parameters(genome_dir, genome_fasta_file, annotations_gff_file):
    # Set up STAR command for genome generation
    star_cmd = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--genomeDir', genome_dir,
        '--genomeFastaFiles', genome_fasta_file,
        '--sjdbGTFfile', annotations_gff_file,
        '--runThreadN', '4',
        '--sjdbGTFtagExonParentTranscript', 'transcript_id',
        '--sjdbGTFfeatureExon', 'exon',
        '--sjdbGTFtagExonParentGene', 'gene_id'
    ]
    
    # Run STAR genome generation
    subprocess.run(star_cmd, check=True)

# Provide the paths for the genome directory, genome FASTA file, and annotations GFF file
genome_dir = '/path/to/genome_dir'
genome_fasta_file = '/path/to/genome.fa'
annotations_gff_file = '/path/to/annotations.gff'

# Generate the genome parameters
generate_genome_parameters(genome_dir, genome_fasta_file, annotations_gff_file)
