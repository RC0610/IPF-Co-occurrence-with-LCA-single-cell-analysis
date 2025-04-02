# Single-Cell RNA Sequencing Analysis: IPF-Co-occurrence-with-LCA/Trajectory and EMT Analysis

## Overview

This repository contains the code and analysis workflow for single-cell RNA sequencing (scRNA-seq) data analysis, focusing on Idiopathic Pulmonary Fibrosis (IPF) and its co-occurrence with lung cancer (LCA). The analysis involves two main parts: **EMT Gene Analysis** and **Trajectory Analysis**. 

### Key Analysis Steps:
1. **EMT Gene Analysis**: Identification of genes related to epithelial-mesenchymal transition (EMT) and saving them for downstream analysis.
2. **Trajectory Analysis**: This part of the analysis involves using Monocle to perform trajectory analysis to infer the developmental paths of different cell types, including pseudotime analysis and differential gene testing based on cell trajectories.

## Directory Structure

- **Trajectory_and_EMT_Analysis**: Contains scripts for EMT gene analysis, trajectory analysis, and result visualization.
- **output**: Contains final output files such as plots, marker genes, differential gene expression results, and processed data.
- **scripts**: Contains R scripts for running the analyses and visualizations.

## Requirements

The following R packages are required for running this analysis:

- Seurat
- tidyverse
- dplyr
- patchwork
- harmony
- data.table
- monocle

To install the required R packages, you can use the following code:

```r
install.packages(c("ggplot2", "patchwork", "dplyr"))
devtools::install_github("satijalab/seurat")
devtools::install_github("immunogenomics/harmony")
devtools::install_github("cole-trapnell-lab/monocle3")
