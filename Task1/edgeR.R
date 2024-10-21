# install.packages("edgeR")

# Load required libraries
library(edgeR)
library(tidyverse)
library(ggplot2)
library(ggrepel)

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
res_tableEG_tb <- topTags(result, n = Inf)$table
res_tableEG_tb <- as.data.frame(res_tableEG_tb)
print(head(res_tableEG_tb))  # Check the first few rows of the results table

# Filter for low expression and significant downregulated genes
edgeR_SM <- res_tableEG_tb %>%
  filter(PValue  < 0.05 & logFC < 0)  # Filtering significant downregulated genes
cat("Number of low expression genes found:", nrow(edgeR_SM), "\n")

# Filter for low expression and significant downregulated genes
edgeR_SM <- res_tableEG_tb %>%
  filter(PValue  < 0.05 & logFC < -0.5)  # Filtering significant downregulated genes
# Print the number of low expression genes found
cat("Number of low expression genes found:", nrow(edgeR_SM), "\n")

# Add a column to indicate significant and low expression genes
res_tableEG_tb <- res_tableEG_tb %>%
  mutate(threshold = case_when(
    PValue < 0.05 & logFC < -0.5 ~ "low_expression_significant",  # Downregulated genes
    PValue < 0.05 & logFC > 0  ~ "high_expression_significant", # Upregulated genes
    TRUE ~ "not_significant"     # Not significant genes
  ))

# Add gene names as a new column in res_tableEG_tb
res_tableEG_tb$gene_name <- rownames(res_tableEG_tb)

# Step 1: Save dispersion estimates plot as PNG
png("Figure/Figure 10. edgeR Volcano plot of SARS gene expression analysis showing differentially expressed genes.png", 
    width = 1600, height = 2000, res = 300)

# Step 2: edgeR volcano plot 
ggplot(res_tableEG_tb, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(colour = threshold), alpha = 0.7) +  # Color points based on significance
  scale_color_manual(values = c("low_expression_significant" = "blue",
                                "high_expression_significant" = "red",
                                "not_significant" = "gray")) +
  geom_text_repel(data = res_tableEG_tb %>%
                    filter(threshold == "low_expression_significant") %>%  # Filter for low expression significant
                    arrange(PValue) %>%                                   # Sort by PValue
                    slice_head(n = 3),                                   # Select top 3
                  aes(label = gene_name),  # Use the new gene_name column
                  max.overlaps = Inf) +  # Allow for all labels
  labs(title = "Volcano Plot of EdgeR Results",
       x = "Log2 Fold Change",
       y = "-Log10 P-Value") +
  theme_minimal() +
  theme(legend.position = "right")

# Step 3: Close the device
dev.off() 

cat("Dataframe with", nrow(edgeR_SM), "rows and", ncol(edgeR_SM), "columns\n")

edgeR_SM
# Optionally, save the results to a CSV file
write.csv(edgeR_SM, "Results_Table/edgeR_SM_low_expression.csv", row.names = TRUE)

#Output Session Information
# Redirect output to a file
sink("Versions_of_all_tools(Session_Information)/sessionInfo_edgeR_R.txt")

# Print sessionInfo() details
sessionInfo()

# Turn off the sink and return output to console
sink()
