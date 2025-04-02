# Cell Interaction Analysis: IPF and LCA

## Overview

This repository contains the code and analysis workflow for investigating cell-cell communication in single-cell RNA sequencing (scRNA-seq) data of Idiopathic Pulmonary Fibrosis (IPF) and its co-occurrence with lung cancer (LCA). The analysis includes:

1. **Data Preparation**: Importing processed data, including cell annotation and expression data, and preparing it for cell-cell interaction analysis.
2. **Cell-Cell Communication**: Using the `CellChat` package to analyze cell communication, focusing on ligand-receptor interactions and signaling pathways.
3. **Network Visualization**: Generating chord diagrams and alluvial plots to visualize the interactions between cell types, with a focus on specific pathways like WNT signaling.
4. **Results Export**: Exporting the interaction matrix and visualizations for downstream analysis and reporting.

## Directory Structure

- **Cell_Interaction_Analysis**: Contains code and analysis for cell-cell communication using the `CellChat` package.
- **output**: Contains the output files, such as interaction results (`interaction_results.txt`), and generated plots (`wnt_in_IPF.pdf` and `interaction_plots_IPF.pdf`).
- **scripts**: Contains the R scripts used for processing, analyzing, and visualizing data.

## Requirements

The following R packages are required for this analysis:

- Seurat
- dplyr
- tidyverse
- data.table
- CellChat
- patchwork
- ggalluvial
- ggsci
- circlize
- RColorBrewer
- ggplot2

To install these packages, run the following R code:

```r
install.packages(c("ggplot2", "patchwork", "dplyr", "tidyverse", "data.table"))
devtools::install_github("jokergoo/CellChat")
devtools::install_github("jokergoo/circlize")
devtools::install_github("hadley/ggalluvial")
