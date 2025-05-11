# JAWesome_RNAproject
This project is part of my NIH-NIDCR K99/R00 Pathway to Independence Award and have several aims. This first READMe contain the scripts and others for completing the first aim of the project, in the form of a manuscript, currently in revision but posted on bioRxv here: https://www.biorxiv.org/content/10.1101/2024.10.02.616385v1.full.

The main objective of the manuscript is to identify and validate tissue-specific (in this case focused on the craniofacial tissues) differentially-expressed genes (DEGs) in a group of dietary specialist pupfishes and overlap these craniofacial-DEGs with genes that we know carry fixed single-nuclear polymorphisms in their regulatory regions. With this pipeline, I analyzed high-depth (>50M reads) tissue-specific (craniofacial vs. caudal tail tissue) RNAseq data from 6 different pupfish species and populations, known for their high craniofacial divergence. We used high-depth RNAseq to identify novel low-abundance transcripts (i.e., chemokine or other membrane receptors) differentially expressed between species and tissues. 

A suggested standard pipeline:
1. Experimental design
2. Isolate the RNA
3. QC the RNA
4. Send to sequence
5. Once you receive the data you will have to:
   a. Check the MF5 files to see if the data was properly downloaded
   b. Do fastqc (QC) check, again
   c. Trimming (if you need to). I used with Trim-galore
   d. Aling your reads to the reference genome with STAR of Kallisto. Kallisto does pseudoalignments and is less sensitive to sequencing depth than STAR. Thus, Kallisto might be more helpful when analyzing low-depth RNAseq data. STAR will also allow you to detect novel splice variants. If you are running this with a genome reference from a model organism, most genome parameters are available from NCBI. If you, like myself, were working with a non-model organism, you will have to generate the genome parameters with the .gff and the genome (.fasta) files with '--runMode genomeGenerate'. The script for this is genomeGenerator_2025.sh
5. Count the aligned transcripts. I used HTSeqCounts()
6. PCA/dendrogram
7. DESeq2
8. Volcano Plot
9. Venn diagram
   
   
