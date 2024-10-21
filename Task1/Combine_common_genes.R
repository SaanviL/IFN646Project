# Load necessary libraries
library(tidyverse)
library(edgeR)
library(DESeq2)
library(tibble)  # For rownames_to_column

# Assuming you have already run DESeq2 and EdgeR analysis and have the result tables

# Load the filtered results from both tools
# Replace these file paths with the actual paths to your results CSV files
deseq2_SM <- read.csv("deseq2_SM_low_expression.csv")
edgeR_SM <- read.csv("edgeR_SM_low_expression.csv", row.names = 1)

# Ensure gene identifiers (rownames) are in the "gene" column for edgeR results
edgeR_SM <- rownames_to_column(edgeR_SM, var = "gene")

# Check the column names in both dataframes to make sure the "gene" column exists
colnames(deseq2_SM)  # Should include "gene"
colnames(edgeR_SM)   # Should include "gene" after rownames_to_column

# Perform the inner join to get common genes between DESeq2 and EdgeR results
combined_results <- inner_join(deseq2_SM, edgeR_SM, by = "gene")
combined_results
# View the number of common genes found and the first few rows of the combined results
cat("Number of common genes found:", nrow(combined_results), "\n")
print(head(combined_results))

# Optionally, save the combined result to a CSV file
write.csv(combined_results, "Results_Table/combined_low_expression_genes.csv", row.names = FALSE)

#Output Session Information
# Redirect output to a file
sink("Versions_of_all_tools(Session_Information)/sessionInfo_Combine_common_genes_R.txt")

# Print sessionInfo() details
sessionInfo()

# Turn off the sink and return output to console
sink()
