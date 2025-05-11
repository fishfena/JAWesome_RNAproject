#!/bin/bash

#SBATCH --job-name=<YOUR.JOB.NAME>
#SBATCH --account=<YOUR.CLUSTER.ACCOUNT>
#SBATCH --time=24:00:00
#SBATCH --partition=<YOUR.PARTITION>
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --array=0-17
#SBATCH --error=<YOUR.ERROR>_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user= <YOUR.USER>


## Do not forget to load STAR as module ##


# Output directory
output_directory="<YOUR.OUTPUT.DIRECTORY>"

# Genome index directory
genome_dir="<GENOME.DIRECTORY>"

# Define input file pairs
read_pairs=(
    "<PATH.TO.YOUR.FILE1>_R1_001.fastq.gz,<PATH.TO.YOUR.FILE1>_R3_001.fastq.gz"
    "<PATH.TO.YOUR.FILE2>_R1_001.fastq.gz,<PATH.TO.YOUR.FILE2>_R3_001.fastq.gz"
	.
	.
	.
    "<PATH.TO.YOUR.FILE_n>_R1_001.fastq.gz,<FILE.TO.YOUR.FILEn>_R3_001.fastq.gz"
)

# Get the current array index
index=$SLURM_ARRAY_TASK_ID

# Extract the paired input files for the current job
input_file1=$(echo "${read_pairs[$index]}" | cut -d',' -f1)
input_file2=$(echo "${read_pairs[$index]}" | cut -d',' -f2)

# Validate input file existence
if [[ ! -f "$input_file1" || ! -f "$input_file2" ]]; then
    echo "ERROR: One or both input files not found:"
    echo "  $input_file1"
    echo "  $input_file2"
    exit 1
fi

# Extract base sample name (this depends on the name of your files)
sample_name=$(basename "$input_file1")
sample_name="${sample_name:0:15}"

# Run STAR
STAR --runThreadN 20 \
--genomeDir "$genome_dir" \
--readFilesIn "$input_file1" "$input_file2" \
--readFilesCommand zcat \
--outFileNamePrefix "${output_directory}/${sample_name}_" \
--outSAMtype BAM SortedByCoordinate \
> "${output_directory}/${sample_name}_STAR.out" \
2> "${output_directory}/${sample_name}_STAR.err"

echo "Alignment completed for ${sample_name}"
