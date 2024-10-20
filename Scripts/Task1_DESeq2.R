# Setup
# Load required libraries
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)

# Load in data
data <- read.table("data/data.txt", header = TRUE, row.names = 1) 
meta <- read.table("meta/meta.txt", header = TRUE, row.names = 1)

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

# Plot PCA
plotPCA(rld, intgroup = "sampletype")


#Hierarchical Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
pheatmap(rld_cor)

## Plot dispersion estimates
dds <- DESeq(dds)
plotDispEsts(dds)


contrast_mk <- c("sampletype", "SARS", "Mock")
levels(dds$sampletype)
res_tableMK_unshrunken <- results(dds, contrast=contrast_mk, alpha = 0.05)

# If needed, install ashr
install.packages('ashr')

res_tableMK <- lfcShrink(dds, contrast=contrast_oe, res=res_tableMK_unshrunken, type="ashr")

plotMA(res_tableMK_unshrunken, ylim=c(-2,2))
plotMA(res_tableMK, ylim=c(-2,2))

class(res_tableMK)

mcols(res_tableMK, use.names=T)

res_tableMK %>% data.frame() %>% View()

### Set thresholds
#baseMean.cutoff<-30


#Filter
res_tableMK_tb <- res_tableMK %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

deseq2_MK <- res_tableMK_tb %>%
  filter(padj  < 0.05 & log2FoldChange < 0)  # Filtering for lower expression levels
cat("Dataframe with", nrow(deseq2_MK), "rows and", ncol(deseq2_MK), "columns\n")
deseq2_MK

#summary(res_tableMK$baseMean)

# Assuming sigMK contains your filtered low expression genes
write.csv(deseq2_MK, "deseq2_MK_low_expression.csv", row.names = FALSE)

