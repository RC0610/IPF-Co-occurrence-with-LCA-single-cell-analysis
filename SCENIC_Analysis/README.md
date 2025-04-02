# SCENIC Analysis for IPF Single-Cell RNA Sequencing Data

## Overview

This repository provides a detailed analysis of single-cell RNA sequencing (scRNA-seq) data from Idiopathic Pulmonary Fibrosis (IPF) using SCENIC (Single-Cell Regulatory Network Inference and Clustering). The analysis focuses on inferring transcriptional regulatory networks, identifying key transcription factors (TFs) and their target genes, and correlating TF activity with gene expression patterns across different cell types.

The main steps of the analysis include:
1. **Data Loading**: Load the single-cell RNA-seq dataset and perform cell type annotation.
2. **Gene Filtering**: Filter genes based on expression levels to retain only informative genes.
3. **Correlation Calculation**: Calculate correlation matrices for gene expression and TF activity.
4. **Genie3 Analysis**: Infer TF-target interactions using the GENIE3 algorithm.
5. **Regulon Inference**: Identify potential regulatory networks (regulons) associated with transcription factors.
6. **AUC Computation**: Compute and binarize AUC scores for each cell and regulon.
7. **Gene and TF Correlation**: Correlate gene expression with TF activity scores.
8. **Visualization**: Generate heatmaps to visualize gene-TF activity correlations.

## Directory Structure

- **SCENIC_Analysis**: Contains the R script for SCENIC analysis, including the steps for gene filtering, correlation calculation, TF-target prediction, and AUC computation.
- **output**: Contains output files including the correlation heatmap (`scenic_plots_IPF_Basal.pdf`), AUC scores, and results from the SCENIC analysis.
- **int**: Intermediate directory for storing results generated during the analysis (e.g., filtered matrices, gene networks).

## Requirements

The following R packages are required for running the SCENIC analysis:

- Seurat
- dplyr
- tidyr
- patchwork
- SCENIC
- harmony
- data.table
- pheatmap
- foreach
- RcisTarget

You can install these packages using the following R code:

```r
install.packages(c("Seurat", "dplyr", "tidyverse", "patchwork", "data.table"))
devtools::install_github("grimpil/SCENIC")
devtools::install_github("satijalab/harmony")
install.packages("pheatmap")
