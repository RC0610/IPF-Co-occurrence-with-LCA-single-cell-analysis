# Single-Cell RNA Sequencing Analysis: IPF-Co-occurrence-with-LCA/Marker_gene_Expression_Analysis_and_Cell_Annotation
## Overview

This repository contains code for the analysis of single-cell RNA sequencing (scRNA-seq) data for IPF (Idiopathic Pulmonary Fibrosis) and its co-occurrence with lung cancer (LCA). The analysis covers multiple steps, including data preprocessing, dimensionality reduction, clustering, differential expression analysis, and cell annotation. 

The workflow begins with reading and processing the pre-processed and clustered data, followed by marker gene identification, cell annotation, and visualization of cell types using UMAP. The key steps involve:
1. **Data prepare**: Loading and preparing the scRNA-seq data that has been pre-processed for analysis.
2. **Cluster Marker Gene Identification**: Identifying marker genes for each cluster using differential expression analysis.
3. **Manual Cell Type Annotation**: Annotating each cluster based on known cell type markers from the literature.
4. **Visualization**: Generating UMAP plots to visualize the annotated cell types in both IPF and LCA datasets.
5. **Results Export**: Exporting the results of differential gene expression and cell annotations.

## Directory Structure

- **Differential_Expression_Analysis_and_Cell_Annotation**: Contains scripts for identifying cluster markers and annotating cell types.
- **output**: Contains output files such as UMAP plots and marker gene lists.
- **scripts**: Contains the R scripts for running the analysis.

## Requirements

To run this analysis, the following R packages are required:

- Seurat
- tidyverse
- dplyr
- patchwork
- harmony
- cowplot
- data.table

You can install these packages using the following R code:

```r
install.packages(c("ggplot2", "patchwork"))
devtools::install_github("satijalab/seurat")
devtools::install_github("immunogenomics/harmony")
