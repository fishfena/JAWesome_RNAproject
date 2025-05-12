# Run this on your personal computer, on the folder you have the .txt counts. This script will merge the counts from line 1 and line 2, for each gene, for each sample. 

# List of data_A and data_B files
data_A_files = ['<YOUR.FILE1_L001>', 
'<YOUR.FILE_L001>',
.
.
.
'<YOUR.FILEn_L001>']

data_B_files = ['<YOUR.FILE2_L001>', 
'<YOUR.FILE_L002>',
.
.
.
'<YOUR.FILEn_L002>']


# Loop through the data files
for data_A_file, data_B_file in zip(data_A_files, data_B_files):
    # Read data from table A (L001)
    with open(data_A_file, 'r') as file_A:
        data_A = [line.strip().split() for line in file_A]

    # Read data from table B (L002)
    with open(data_B_file, 'r') as file_B:
        data_B = [line.strip().split() for line in file_B]

    # Combine data from both tables and Calculate the sum
    total_data = []
    for row_A, row_B in zip(data_A, data_B):
        gene_A, count_reads_A = row_A
        gene_B, count_reads_B = row_B

        # Verify that the genes match between the two tables
        if gene_A != gene_B:
            raise ValueError("Gene names do not match between tables A and B")

        count_reads_A = int(count_reads_A)
        count_reads_B = int(count_reads_B)
        total_reads = count_reads_A + count_reads_B

        total_data.append([gene_A, count_reads_A, gene_B, count_reads_B, total_reads])

    # Extract the first 6 characters of the data_A filename
    file_prefix = data_A_file[:6]

    # Create the output filename
    output_filename = f'{file_prefix}_total_reads.txt'

    # Write the combined data to the new table
    with open(output_filename, 'w') as total_file:
        total_file.write("gene\tcount_reads_L001\tgene\tcount_reads_L002\tcount_reads_L001+L002\n")
        for row in total_data:
            total_file.write("\t".join(map(str, row)) + "\n")

    print(f"New table '{output_filename}' has been created.")
