# Single-Cell RNA Sequencing Analysis: IPF-Co-occurrence-with-LCA/clustering analysis

## Overview

This repository contains the code and analysis workflow for single-cell RNA sequencing (scRNA-seq) analysis of IPF (Idiopathic Pulmonary Fibrosis) and its co-occurrence with lung cancer (LCA). The analysis focuses on data preprocessing, quality control (QC), dimensionality reduction, clustering, and visualization of the results.

The primary analysis includes:

1. **Data prepares**: Reads data that has been applied quality control filters, and prepares the data for downstream analysis.
2. **Dimensionality Reduction and Clustering**: Uses PCA and Harmony batch correction, followed by clustering of cells.
3. **Visualization**: Visualizes the results using UMAP (Uniform Manifold Approximation and Projection).
4. **Results Export**: Outputs clustering results and other related files for further analysis.

## Directory Structure

- **Clustering_Analysis**: Contains code for dimensionality reduction, clustering, and results visualization.
- **output**: Contains final output files, such as UMAP plots, clustering results, and processed data.

## Requirements

The following R packages are required for analysis:

- Seurat
- ggplot2
- patchwork
- harmony

You can install these packages using the following R code:

```r
install.packages(c("ggplot2", "patchwork"))
devtools::install_github("satijalab/seurat")
devtools::install_github("immunogenomics/harmony")
