# JAWesome_RNAproject
I will use this pipeline to analyze the high-depth RNAseq data from the JAWesome project. 
The first step of this project was to isolate RNA from the tissue. The tissue that we wanted to compare for differential expression was the branchial apparatus (branchial arches) and the jaw (lower and upper jaw) with the tail (a tissue that hasn't diverged phenotypically) between different species of pupfishes endemic to San Salvador Island, Bahamas, and other pupfishes as outgroup species. 
This pipeline will be highly detailed, allowing anyone to do RNAseq analyses from scratch. There were a lot of things that if have been out there in simple terms, would have increased the processing time by A LOT. 
Here's the pipeline. Scripts are available in the Script folder. 

1. Isolate the RNA
2. QC the RNA
3. Send to your preferred sequencing facility
4. Once you receive the data you will have to:
   a. Check the MF5 files for checking that the data was properly downloaded
   b. Do fastqc (QC) check
   c. Trimming (if you need to) with Trim-galore
   d. Aling your reads to the reference genome with STAR of Kallisto. Kallisto does pseudoalignment and is less sensitive to sequencing depth than STAR. Thus, Kallisto might be more helpful when analyzing    low-depth RNAseq data. STAR will also allow you to detect novel splice variants. If you are running this with a genome reference from a model organism, most genome parameters are available from NCBI.      If you, like myself, were working with a non-model organism, you will have to generate the genome parameters with the .gff and the genome (.fasta) file with '--runMode genomeGenerate'. The script for      this is genome_dir_script.py
   e. Next: featurecounts()
   
