# JAWesome_RNAproject
I will use this pipeline to analyze the high-depth RNAseq data from the JAWesome project. The main objectives of this project are:
1. Compare the RNA profile from the jaw between species and identify DE genes
2. Compare the RNA profile from the jaw vs tail across pupfishes and identify DE genes
3. Identify DE low and high-abundance transcripts from (1) and (2)
4. Identify DE splice variants from (1) and (2)
5. Build a network of interacting genes from (1) and (2)

The first step of this project was to isolate RNA from the tissue. The tissue that we will compare for DE was the branchial apparatus (branchial arches) and the jaw (lower and upper jaw) with the tail (a tissue that hasn't diverged phenotypically) between different species of pupfishes endemic to San Salvador Island, Bahamas, and other pupfishes as outgroup species. Biological replicates: minimum of 3 and maximum of 5, per species. Each biological replicate consists of the dissected tissues from 5-7 larvae. 

Scripts are available in the Script folder.

Here's the pipeline up to date, Oct 8th, 2024.

1. Isolate the RNA. DM-me if details are needed
2. QC the RNA, in this case, w/ Bioanalyzer
3. Send to your preferred sequencing facility
4. Once you receive the data you will have to:
   a. Check the MF5 files for checking that the data was properly downloaded
   b. Do fastqc (QC) check
   c. Trimming (if you need to) with Trim-galore
   d. Aling your reads to the reference genome with STAR of Kallisto. Kallisto does pseudoalignments and is less sensitive to sequencing depth than STAR. Thus, Kallisto might be more helpful when analyzin  	low-depth RNAseq data. STAR will also allow you to detect novel splice variants. If you are running this with a genome reference from a model organism, most genome parameters are available from NCBI.      If you, like myself, were working with a non-model organism, you will have to generate the genome parameters with the .gff and the genome (.fasta) files with '--runMode genomeGenerate'. The script for      this is genome_dir_script.py
   e. HTSeqCounts()
   f. PCA/dendrogram
   g. DESeq2
   h. Volcano Plot
   i. Venn diagram
   
   
