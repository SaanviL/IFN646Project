# IFN646 - COVID-19 Gene Expression and CRISPR Guide Design
(Group 9)

# Introduction
This project aims to identify genes with significantly lower expression in COVID-19 patients and design suitable CRISPR guide RNA (gRNA) sequences to potentially restore the deregulated gene expression. The project involves multiple stages, from analyzing RNA-seq data to evaluating CRISPR gRNA designs across different human genomes.

This README provides an overview of the data and instructions on how to run the analysis.

# Data
The data required for this project includes RNA-seq raw counts, human genome sequences, annotations, and variant information. Each part includes data and outputs corresponding to the specific file names listed below.

**Datasets used can be found in the 'Data' folder in Task 1 and Task_3 Data folder.**

## Task 1: Gene Expression Analysis Data
Task 1 focuses on analyzing RNA-seq raw counts to identify genes with significantly lower expression in COVID-19 patients using the DESeq2 and edgeR packages. All data, code, and output for Task 1 are located in the **`Task1`** folder.

### **Task1** Folder Included:
### Folder:
- Data: **`data.txt"`** - RNA-seq raw counts data.
- Figure: Output figures generated from the DESeq2 and edgeR scripts.
- Meta: **`meta.txt"`** - Information regarding sample types used in the analysis.
- Results Table: Output table of significant low-expression genes generated from the following R scripts:
1. **`DESeq2.R`**: **`deseq2_SM_low_expression.csv`**
2. **`edgeR.R`**: **`edgeR_SM_low_expression.csv`**
3. **`Combine_common_genes.R`**: **`combined_low_expression_genes.csv`**
- Version Information: Session information for the three R scripts to ensure reproducibility.
1. **`DESeq2.R`**: **`sessionInfo_DESeq2_R.txt`**
2. **`edgeR.R`**: **`sessionInfo_edgeR_R.txt`**
3. **`Combine_common_genes.R`**: **`sessionInfo_Combine_common_genes_R.txt`**

### Code:
1. RStudio Project: Task1.Rproj - Organizes all related scripts and resources.
2. R Script: **`DESeq2.R`** - Code for finding low-expression genes using the DESeq2 package.
3. R Script: **`edgeR.R`** - Code for finding low-expression genes using the edgeR package.
4. R Script: **`Combine_common_genes.R`** - Code for combining DESeq2 and edgeR results and selecting intersecting genes as the final list of significant low-expression genes.


## Task 3: Variant Matching and gRNA Design
Task 3 involves matching variants from individuals on chromosome 1 with variants on chromosome 2 within the off-target position range to identify potential gRNA targets. All data, code, and output for Task 3 are located in the **`Task3_Matched_Code`** folder.

### **Task3_Matched_Code** Folder Included:
- Task_3_Data:
1. Off-Target Table: **`off_target_table.csv`** - List of 38/rev guide off-targets.
2. Variants from Chromosome 1: **`sampled_chr1.vcf`** - Variant data from 7 individuals on chromosome 1.
3. Variants from Chromosome 2: **`sampled_chr2.vcf`** - Variant data from the same individuals on chromosome 2.
- Matched Genes Table: **`matched_genes.csv`** Variant data of matched genes.

#### Code:
- Python Code: **`IFN_464_Project_Task3.py`** - Code for matching variants and analyzing off-targets.
- Jupyter Notebook: **`IFN_464_Project_Task3.ipynb`** - Jupyter version of the same code for interactive analysis.

#Other:
- README for Task 3:Instructions for using **`IFN_464_Project_Task3.py and **`IFN_464_Project_Task3.ipynb`**.

# Libraries
DESeq2 and edgeR

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
```


