#!/bin/bash

#SBATCH --job-name=<YOUR.JOB.NAME>
#SBATCH --account=<YOUR.CLUSTER.ACCOUNT>
#SBATCH --time=24:00:00
#SBATCH --partition=<YOUR.PARTITION>
#SBATCH --error=<YOUR.ERROR>_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user= <YOUR.USER>

# Create output directory for STAR genome index
mkdir -p <GENOME.DIRECTORY>

# Run STAR genome generation
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir <GENOME.DIRECTORY> \
--genomeFastaFiles <PATH.TO.FASTA.GENOME.DIRECTORY> \
--sjdbGTFfile <PATH.TO.GTF.GENOME.DIRECTORY> \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript 'transcript_id' \
--sjdbGTFfeatureExon 'exon' \
--sjdbGTFtagExonParentGene 'gene_id'
