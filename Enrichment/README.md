## Single-Cell RNA Sequencing Analysis: IPF-Co-occurrence-with-LCA/Enrichment Analysis
## Overview

This analysis pipeline performs Gene Ontology (GO) Biological Process (BP) and KEGG pathway enrichment analysis using a set of differentially expressed genes (DEGs) derived from various sources. The analysis includes the extraction of genes from multiple datasets related to IPF and its co-occurrence with lung cancer (LCA), followed by GO and KEGG enrichment analysis. The goal of this analysis is to identify biological processes and pathways associated with IPF and its different cell types.

The primary analysis steps include:
1. **Gene Extraction**: Extract genes from different datasets such as time-difference analysis, BEAM analysis, monocle DEGs, EMT gene set, and IPF filtered DEGs.
2. **Intersection of Genes**: Identify the common genes across all datasets to focus on the most relevant genes for enrichment analysis.
3. **GO BP Enrichment Analysis**: Perform GO Biological Process enrichment analysis on the common genes.
4. **KEGG Pathway Enrichment Analysis**: Perform KEGG pathway enrichment analysis using the common genes.
5. **Visualization**: Generate bubble plots to visualize the top enriched pathways and biological processes, and save the results for downstream analysis.

## Directory Structure

- **Enrichment**: Contains the code for GO and KEGG enrichment analysis, as well as the related gene files.
- **output**: Contains the results of enrichment analysis, including bubble plots for GO and KEGG pathways and text files with detailed enrichment results.

## Requirements

To run this analysis, the following R packages are required:

- **dplyr**
- **tidyverse**
- **DESeq2**
- **msigdbr**
- **GSEABase**
- **clusterProfiler**
- **org.Hs.eg.db**
- **enrichplot**
- **data.table**
- **ggplot2**

These packages can be installed using the following commands in R:

```r
install.packages(c("dplyr", "tidyverse", "DESeq2", "GSEABase", "clusterProfiler", "enrichplot", "data.table", "ggplot2"))
devtools::install_github("msigdbr/msigdbr")
devtools::install_github("immunogenomics/harmony")
