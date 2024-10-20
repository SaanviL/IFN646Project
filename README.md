# IFN646 - COVID-19 Gene Expression and CRISPR Guide Design Project (Group 9)

# Introduction
This project aims to identify genes with significantly lower expression in COVID-19 patients and design suitable CRISPR guide RNA (gRNA) sequences to potentially restore the deregulated gene expression. The project involves multiple stages, from analyzing RNA-seq data to evaluating CRISPR gRNA designs across different human genomes.

This README provides an overview of the data and instructions on how to run the analysis.

# Data
The data required for this project includes RNA-seq raw counts, human genome sequences, annotations, and variant information.

**Datasets used can be found in the 'Data' folder. Larger data files used is available with the following link.**

Human Genome assembly GRCh38 = https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/

**- Genomic data (GCF_000001405.26_GRCh38_genomic.fna.gz, 902MB compressed)**

**- Genome annotation (GCF_000001405.26_GRCh38_genomic.gff.gz, 24MB compressed)**

# Libraries
DESeq2 and edgeR

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
```


