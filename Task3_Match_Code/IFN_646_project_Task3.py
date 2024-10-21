# -*- coding: utf-8 -*-
"""IFN_646_Project_Task3.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/17w91T_qBAbqLMavhNnS0agmkv56DrI9R
"""

import pandas as pd

# Load the off-target data
off_target_data = pd.read_excel("Task3_Data/off_target_table.xls", header=None)
off_target_data.columns = off_target_data.iloc[8]  # Set column names from the 8th row
off_target_data = off_target_data.iloc[9:]  # Keep the data rows

# Optional: Reset the index after filtering
off_target_data.reset_index(drop=True, inplace=True)

# Filter for guideId: '1rev' and chromosomes 'chr1' or 'chr2'
filtered_off_target_data = off_target_data[(off_target_data['guideId'] == '1rev') &
                                           (off_target_data['chrom'].isin(['chr1', 'chr2']))]

# Define a mapping for VCF files based on chromosome
vcf_file_mapping = {
    'chr1': 'Task3_Data/sampled_chr1.vcf',
    'chr2': 'Task3_Data/sampled_chr2.vcf',
}

# Initialize a list to store matched genes
matched_genes = []

# Iterate through each row in the filtered off-target data
for index, row in filtered_off_target_data.iterrows():
    chrom = row['chrom'].strip()  # This will be 'chr1' or 'chr2'
    start = int(row['start'])
    end = int(row['end'])

    # Strip 'chr' from chrom to match the VCF #CHROM format
    chrom_num = chrom.replace('chr', '').strip()  # Convert 'chr1' to '1', 'chr2' to '2'

    # Determine which VCF file to use based on the chromosome
    vcf_file = vcf_file_mapping.get(chrom)  # Get the VCF file based on chromosome

    if vcf_file:  # Proceed only if a valid VCF file is found
        # Load the VCF file
        vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)  # Modify according to VCF format
        vcf_data.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                            'HG00142', 'HG00178', 'HG00237', 'HG01521', 'HG03129', 'NA19175', 'NA19462']

        # Ensure 'POS' and '#CHROM' are correctly formatted
        vcf_data['POS'] = vcf_data['POS'].astype(int)  # Convert 'POS' to integers
        vcf_data['#CHROM'] = vcf_data['#CHROM'].astype(str).str.strip()  # Remove any extra spaces from '#CHROM'


        # Filter the VCF data for positions within the start and end range, and matching chromosome
        filtered_vcf = vcf_data[(vcf_data['#CHROM'] == chrom_num) & (vcf_data['POS'] >= start) & (vcf_data['POS'] <= end)]


        # Print matching rows in the VCF file
        if not filtered_vcf.empty:
            print(f"Matches found for guideId {row['guideId']} in {vcf_file}:")
            matched_genes.append(filtered_vcf)

            print(matched_genes)

            print("\n")  # New line for better readability

        else:
            print(f"No matches found for guideId {row['guideId']} in {vcf_file} within range {start}-{end}.")
    else:
        print(f"No VCF file found for chromosome {chrom}.")

# Combine all matched VCF DataFrames into one DataFrame
if matched_genes:  # Check if there are any matched genes
    matched_genes_df = pd.concat(matched_genes, ignore_index=True)

    # Save the matched_genes DataFrame to a CSV file
    output_file = "matched_genes.csv"
    matched_genes_df.to_csv(output_file, index=False)
    print(f"There are {len(matched_genes_df)} matched genes have been saved to {output_file}.")  # Confirmation message
else:
    print("No matched genes found.")

matched_genes_df

