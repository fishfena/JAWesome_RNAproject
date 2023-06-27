# JAWesome_RNAproject
I will use this pipeline to analyze the high-depth RNAseq data from the JAWesome project. 
The first step of this project was to isolate RNA from the tissue. The tissue that we wanted to compare for differential expression was the branchial apparatus (branchial arches) and the jaw (lower and upper jaw) with the tail (a tissue that hasn't diverged phenotypically) between different species of pupfishes endemic to San Salvador Island, Bahamas, and other pupfishes as outgroup species. 
This pipeline will be highly detailed, allowing anyone to do RNAseq analyses from scratch. There were a lot of things that if have been out there in simple terms, would have increased the processing time by A LOT. 

1. Isolate the RNA
2. QC the RNA
3. Send to your preferred sequencing facility
4. Once you receive the data you will have to:
   a. Check the MF5 files for checking that the data was properly downloaded
   b. Do fastqc (QC) check
   c. Trimming (if you need to) with Trim-galore
   d. Aling your reads to the reference genome
   
