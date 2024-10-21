# Setup
# Load required libraries
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(ggplot2)
library(ggrepel)

# Load in data
data <- read.table("Data/data.txt", header = TRUE, row.names = 1) 
meta <- read.table("Meta/meta.txt", header = TRUE, row.names = 1)

# Check classes of the data we just brought in
class(data)
class(meta)

# Check that sample names match in both files
all(colnames(data) %in% rownames(meta)) # TRUE if all samples in data are in meta
all(colnames(data) == rownames(meta))    # TRUE if sample names match

# Step 1: Subset the count data for SARS and Mock samples
count_data_subset <- data[, c("S_2_SARS_5dpi_S70003", "S_3_SARS_5dpi_S69996", "S_2_mock_5dpi_S70002", "S_3_mock_5dpi_S69997")]

# Step 2: Define the condition factor
sampletype <- factor(c("SARS", "SARS", "Mock", "Mock"))

# Step 3: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data_subset, 
                              colData = data.frame(sampletype), 
                              design = ~ sampletype)

# Quality Control
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Normalization
normalized_counts <- counts(dds, normalized = TRUE)

# Transform counts for data visualization
rld <- rlog(dds, blind = TRUE)

#Output Plot
# Step1: Save PCA plot as PNG
png("Figure/Figure 1. Principal components analysis (PCA) of SARS and Mock Samples.png", 
    width = 1600, height = 1200, res = 300)

# Step2: Plot PCA 
plotPCA(rld, intgroup = "sampletype")

# Step 3: Close the device
dev.off()

#Hierarchical Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor) 

#Output Plot
# Step1: Save Heatmap plot as PNG
png("Figure/Figure 2. Hierarchical Clustering Heatmap of SARS and Mock Samples.png", 
    width = 1600, height = 1200, res = 300)

# Step2:Plot heatmap
pheatmap(rld_cor)

# Step 3: Close the device
dev.off()

dds <- DESeq(dds)

# Step 1: Save dispersion estimates plot as PNG
png("Figure/Figure 3. DESeq2 Dispersion Estimate for SARS.png", 
    width = 1600, height = 1200, res = 300)

# Plot dispersion estimates
plotDispEsts(dds)  # This generates the plot

# Step 3: Close the device
dev.off()  # This saves the plot to the specified file

# contrast
contrast_SM <- c("sampletype", "SARS", "Mock")
levels(dds$sampletype)
res_tableDES_unshrunken <- results(dds, contrast=contrast_SM, alpha = 0.05)

# If needed, install ashr
install.packages('ashr')

res_tableDES <- lfcShrink(dds, contrast=contrast_sm, res=res_tableDES_unshrunken, type="ashr")

# Step 1: Save dispersion estimates plot as PNG
png("Figure/Figure 4. MA Plots for Differential Expression Analysis. (A) Unshrunken estimates.png", 
    width = 1600, height = 2000, res = 300)

# Step2: plotunshrunken estimates
plotMA(res_tableDES_unshrunken, ylim=c(-2,2))

# Step 3: Close the device
dev.off() 

# Step 1: Save dispersion estimates plot as PNG
png("Figure/Figure 4. MA Plots for Differential Expression Analysis. (B) Shrunken estimates.png", 
    width = 1600, height = 2000, res = 300)

# Step2: plot shrunken estimates
plotMA(res_tableDES, ylim=c(-2,2))

# Step 3: Close the device
dev.off() 

class(res_tableDES)

mcols(res_tableDES, use.names=T)

res_tableDES %>% data.frame() %>% View()

res_tableDES_tb <- res_tableDES %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#Volcano Plot


# IF we set optimal number pvlause<0.05 and log2fold change in < 0 => only 2 gene meet this requirement
plot_data <- res_tableDES_tb %>%
  mutate(threshold_SM = case_when(
    padj > 0.05 & log2FoldChange > 0 ~ "high_expression_significant",  # Red
    padj < 0.05 & log2FoldChange < 0 ~ "low_expression_significant",   # Another color, e.g., blue
    TRUE ~ "not_significant"
  ))

# Step1: Output the plot to a PNG file
png("Figure/ Figure 5.DESeq2 Volcano plot of SARS gene expression with default thresholds showing differentially expressed genes.png", 
    width = 1600, height = 2000, res = 300)

# Step2: Volcano plot with color mapping
ggplot(plot_data) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_SM)) +
  ggtitle("SARS Expression Analysis") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_color_manual(values = c(
    high_expression_significant = "red",      # Assign red for high expression significant
    low_expression_significant = "blue",      # Assign blue for low expression significant
    not_significant = "gray"                  # Assign gray for not significant
  )) +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Step 3: Close the device
dev.off() 

# Threshold Experiment  
deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.05 & log2FoldChange < 0)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# => Dataframe with 2 rows and 7 columns

# Therefore, So we relax the p-value to see to generate more gene  
deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.1 & log2FoldChange < 0)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# => Dataframe with 7 rows and 7 columns

deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.15 & log2FoldChange < 0)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# => Dataframe with 15 rows and 7 columns

deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.2 & log2FoldChange < 0)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# => Dataframe with 24 rows and 7 columns 

# P-value over 0.2 can have 20 genes so we set p-value<0.2 and try to find more less expression level so to change tight log2 fold

deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.2 & log2FoldChange < -0.1)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# =>- Dataframe with 22 rows and 7 columns

deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.2 & log2FoldChange < -0.2)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# =>- Dataframe with 22 rows and 7 columns

deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.2 & log2FoldChange < -0.5)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# => Dataframe with 18 rows and 7 columns

#Finally: 
deseq2_SM <- res_tableDES_tb %>%
  filter(padj  < 0.2 & log2FoldChange < -0.5)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_SM), "rows and", ncol(deseq2_SM), "columns\n")
# => Dataframe with 18 rows and 7 columns
deseq2_SM

# Volcano

# Step1: Output the plot to a PNG file
png("Figure/ Figure 7.DESeq2 Volcano plot of SARS gene expression with the modified adjusted p-value and log2FoldChange thresholds showing differentially expressed genes.png", 
    width = 1600, height = 2000, res = 300)

# Create a separate dataframe for plotting with thresholds
plot_data <- res_tableDES_tb %>%
  mutate(threshold_SM = case_when(
    padj > 0.05 & log2FoldChange > 0 ~ "high_expression_significant",  # Red
    padj < 0.2 & log2FoldChange < -0.5 ~ "low_expression_significant",   # Another color, e.g., blue
    padj > 0.05 ~ "not_significant"
  ))

# Volcano plot with color mapping and labels for the top 3 low expression significant genes
ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_SM), alpha = 0.7) +  # Color points based on significance
  geom_text_repel(data = plot_data %>% 
                    filter(threshold_SM == "low_expression_significant") %>%  # Filter for low expression significant
                    arrange(padj) %>%                                          # Sort by adjusted p-value
                    slice_head(n = 3),                                       # Select top 3
                  aes(label = gene),  # Use the correct label column
                  na.rm = TRUE,  # Remove NA values
                  max.overlaps = Inf) +  # Allow for all labels
  ggtitle("SARS Expression Analysis") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted P-Value") +
  scale_color_manual(values = c(
    high_expression_significant = "red",      # Assign red for high expression significant
    low_expression_significant = "blue",      # Assign blue for low expression significant
    not_significant = "gray"                  # Assign gray for not significant
  )) +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

# Assuming sigMK contains your filtered low expression genes
write.csv(deseq2_SM, "Results_Table/deseq2_SM_low_expression.csv", row.names = FALSE)

#Output Session Information
# Redirect output to a file
sink("Versions_of_all_tools(Session_Information)/sessionInfo_DESeq2_R.txt")

# Print sessionInfo() details
sessionInfo()

# Turn off the sink and return output to console
sink()
