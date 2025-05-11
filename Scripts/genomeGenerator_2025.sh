#!/bin/bash

#SBATCH --job-name=generateGenome_variegatus_nonabitio
#SBATCH --account=fc_fishes
#SBATCH --time=24:00:00
#SBATCH --partition=savio3
#SBATCH --error=generateGenome_nonabitio3.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=mfpalominos@berkeley.edu

# Create output directory for STAR genome index
mkdir -p cVariegatusSTAR_nonabitio3

# Run STAR genome generation
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir cVariegatusSTAR_nonabitio3 \
--genomeFastaFiles /global/scratch/users/mfpalominos/RNAseqdata/genomeDir/Cyprinodon_variegatus.C_variegatus-1.0.dna.toplevel.fa \
--sjdbGTFfile /global/scratch/users/mfpalominos/RNAseqdata/genomeDir/Cyprinodon_variegatus.C_variegatus-1.0.110.gtf \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript 'transcript_id' \
--sjdbGTFfeatureExon 'exon' \
--sjdbGTFtagExonParentGene 'gene_id'
