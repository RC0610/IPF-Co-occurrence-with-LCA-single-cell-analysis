Single-Cell RNA-seq Analysis of IPF and LCA Co-occurrence
This project aims to explore the co-occurrence of Idiopathic Pulmonary Fibrosis (IPF) and Lung Cancer (LCA) through single-cell RNA sequencing (scRNA-seq) analysis, with a particular focus on Basal Cells and Type 2 Alveolar Epithelial Cells in IPF and their relationship to LCA. The analysis covers a range of methodologies including data preprocessing, quality control, clustering analysis, marker gene analysis, differential gene expression, trajectory analysis, EMT gene analysis, functional enrichment analysis, cell interaction analysis, and copy number variation (CNV) analysis.

Project Structure
The project is organized into several modules, each corresponding to a specific step in the analysis pipeline.

/Preprocessing_and_QC: This folder contains scripts for data preprocessing and quality control (QC). It includes a quality control report (QC_report.txt) and visualizations (QC_plot.pdf).

/Clustering_Analysis: Contains scripts for clustering analysis and dimensionality reduction. Results are stored in clustering_results.txt and visualized in clustering_plot.pdf.

/Marker_gene_Expression_Analysis_and_Cell_Annotation: This folder focuses on identifying marker genes and annotating cell types. The results are stored in filtered_marker_genes.txt and visualized in annotation_plot.pdf.

/Differential_Expression_Analysis: This module performs differential gene expression analysis, with results saved in filtered.Degs.txt.

/Trajectory_and_EMT_Analysis: Focuses on analyzing cell trajectories and EMT-related genes. Results are saved in trajectory_results.txt, trajectory_plot.pdf, and EMT_gene_results.txt.

/Enrichment: This folder performs functional enrichment analysis, and results are stored in enrichment_results.txt.

/Cell_Interaction_Analysis: Contains scripts for cell interaction analysis and functional enrichment, with results stored in interaction_results.txt and visualized in interaction_plots.pdf.

/SCENIC_Analysis: This module performs SCENIC analysis for gene regulatory network inference. Results are stored in scenic_plots_IPF_Basal.pdf and intermediate files are located in the scenic_results/ folder.

/CNV_Analysis: Includes CNV analysis with results stored in cnv_results.txt and visualizations in cnv_plots.pdf.

/Data: Contains raw and processed data files such as scRNA.RDS (raw data) and scRNA_harmony_clustered.RDS (processed data).

Running the Analysis
To run the analysis, clone the repository and ensure that the necessary R packages are installed. The main dependencies for the analysis are:

dplyr

tidyverse

patchwork

Seurat

infercnv

SCENIC

Each module contains its own script, which can be run independently. The results for each analysis step are stored in text files or PDF plots, making it easy to review and visualize the analysis outcomes.