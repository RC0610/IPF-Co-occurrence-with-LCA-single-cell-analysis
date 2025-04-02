# CNV Analysis for LCA Single-Cell RNA Sequencing Data

This section describes the workflow for performing **Copy Number Variation (CNV)** analysis on single-cell RNA sequencing (scRNA-seq) data of **Lung Cancer (LCA)**. The analysis uses the `infercnv` R package to infer CNV patterns in single-cell RNA-seq data.

## Prerequisites

- **R packages**: 
  - `Seurat`
  - `dplyr`
  - `tidyverse`
  - `infercnv`
  - `patchwork`
  - `data.table`
  - `ggpubr`
  - `RColorBrewer`

- **Input Files**:
  - `scRNA_LCA_harmony_clustered.RDS`: Seurat object containing the clustered scRNA-seq data for LCA.
  - `ann_CA.csv`: Cell annotation file for LCA samples.
  - `human.gene.positions`: File containing human gene position information used for CNV analysis.
  - `groupFiles.txt`: A tab-delimited text file containing the groupings of the cells (e.g., cell types).
  
- **Output Files**:
  - `geneLocate.txt`: Processed file with gene positions.
  - `groupFiles.txt`: File with cell type annotations for CNV analysis.
  - `cnv boxplot.pdf`: A boxplot of CNV scores across different cell types.
  - `infercnv.pdf`: CNV heatmap visualizing CNV across cells.

## Workflow

### 1. **Load Data**:
```r
scRNA = readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_LCA_harmony_clustered.RDS")
celltype = fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation/ann_CA.csv")
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)) {
    scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]
}
