#Single-Cell RNA Sequencing Analysis: IPF-Co-occurrence-with-LCA/Differential_Expression_Analysis
**Analysis Workflow**  
This code performs differential gene expression analysis on single-cell RNA sequencing data for Idiopathic Pulmonary Fibrosis (IPF) and Lung Cancer Association (LCA) cohorts. Key steps include:  
1. Environment setup and package loading  
2. Loading pre-processed Seurat objects clustered via Harmony(scRNA_IPF_harmony_clustered.RDS/scRNA_LCA_harmony_clustered.RDS)  
3. Cell type annotation using marker genes from literature (ann_IPF.csv/ann_CA.csv)  
4. Wilcoxon rank-sum test for cell-cluster marker identification  
5. Filtering criteria: adj.p < 0.05, log2FC > 0.25, pct.1 > 25%, pct.2 < 10%  
6. Export of significant markers to filtered.Degs.IPF/LCA.txt  

**Requirements**  
```r
# Install dependencies
install.packages(c("Seurat", "tidyverse", "dplyr", "patchwork", "cowplot", "data.table", "ggplot2"))
BiocManager::install("SingleR")
remotes::install_github("immunogenomics/harmony")