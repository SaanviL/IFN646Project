# install.packages("edgeR")

# Load required libraries
library(edgeR)
library(tidyverse)

# Load in data
data <- read.table("data/data.txt", header = TRUE, row.names = 1) 
meta <- read.table("meta/meta.txt", header = TRUE, row.names = 1)

# Subset the count data for the samples of interest
count_data_subset <- data[, c("S_2_SARS_5dpi_S70003", "S_3_SARS_5dpi_S69996", "S_2_mock_5dpi_S70002", "S_3_mock_5dpi_S69997")]

# Define sample types
sampletype <- factor(c("SARS", "SARS", "Mock", "Mock"))

# Create DGEList object
dge <- DGEList(counts = count_data_subset, group = sampletype)

# Normalize data
dge <- calcNormFactors(dge)

# Define the design matrix
design <- model.matrix(~ sampletype)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmFit(dge, design)

# Perform likelihood ratio test
result <- glmLRT(fit)  # Change coef depending on your comparison of interest

# Extract results
res_table <- topTags(result, n = Inf)$table
res_table <- as.data.frame(res_table)
print(head(res_table))  # Check the first few rows of the results table

# Filter for low expression and significant downregulated genes
edgeR_MK <- res_table %>%
  filter(PValue  < 0.05 & logFC < 0)  # Filtering significant downregulated genes

# Print the number of low expression genes found
cat("Number of low expression genes found:", nrow(edgeR_MK), "\n")
edgeR_MK
# Optionally, save the results to a CSV file
write.csv(edgeR_MK, "edge_MK_low_expression.csv", row.names = TRUE)

