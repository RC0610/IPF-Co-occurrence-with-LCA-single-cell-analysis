# Preprocessing and QC for IPF and LCA Single-Cell RNA-Seq Data/Preprocessing_and_QC

This directory contains the scripts and functions used for preprocessing and quality control (QC) of single-cell RNA sequencing (scRNA-seq) data from IPF (Idiopathic Pulmonary Fibrosis) and LCA (Lung Cancer with IPF co-occurrence). Below, we describe the various steps involved in the data processing pipeline.

## Directory Structure

- **data**: Contains the raw data for IPF datasets.
  - **h5**: Directory containing the `.h5` files for the GSE143706 dataset.
  - **matrix1**: Directory containing matrix files for the GSE128033 dataset.
  - **matrix2**: Directory containing matrix files for the GSE190889 dataset.
- **data**: Contains the raw data for LCA datasets.
  - **matrix1**: Directory containing matrix files for the GSE198099 dataset.
  - **matrix2**: Directory containing matrix files for the GSE162498 dataset.
 
## Steps in Data Preprocessing and Quality Control

### 1. **Data Loading**
The first step involves loading the raw data from different sources:
- **GSE143706**: Loaded from `.h5` files.
- **GSE128033** and **GSE190889**: Loaded from matrix directories.
- **GSE198099** and **GSE162498**: Loaded from matrix directories.

The `read_sc_h5` and `read_sc_matrix` functions are used to load the respective datasets, with each sample assigned a unique sample ID.

### 2. **Data Merging**
After loading the individual datasets, the data from each source (IPF or LCA) are merged into a single Seurat object for further processing.

### 3. **Quality Control (QC)**
To ensure the data quality, several QC steps are performed:
- **Gene count per cell (nFeature_RNA)**: Cells with too few or too many genes are filtered out.
- **UMI count per cell (nCount_RNA)**: Cells with unusually low or high UMI counts are removed.
- **Mitochondrial content (percent.mt)**: Cells with high mitochondrial gene content are excluded, as this may indicate dying or damaged cells.
- **Hemoglobin content (percent.HB)**: Cells with high hemoglobin gene content are removed to avoid contamination from erythrocytes.

The QC thresholds used are based on recommendations from Satija Lab and other studies. The filtering process removes low-quality cells and ensures that only high-quality cells are retained for further analysis.

### 4. **Data Filtering**
The data is filtered using the criteria mentioned above. Cells that do not pass the thresholds are removed, and the data is saved before and after filtering. QC plots are generated to visualize the quality control metrics.

### 5. **Normalization**
The data is normalized to account for differences in sequencing depth across cells. Seurat’s normalization function is used to scale the data, ensuring that expression levels are comparable across cells.

### 6. **Harmonization**
To remove any batch effects between the different datasets, Harmony is applied for batch correction. This ensures that any dataset-specific variation does not obscure biological variation.

### 7. **Cell Cycle Scoring**
Cell cycle phase scoring is performed to regress out the effects of cell cycle variation on downstream analyses.

### 8. **Final Data Saving**
After the preprocessing and QC steps, the final processed data is saved in an RDS file (`scRNA_IPF_QC_Harmony_processed.RDS`). If LCA data is being processed, the file name changes to `scRNA_LCA_QC_Harmony_processed.RDS`.

## QC Report
A detailed QC report is generated for each dataset, summarizing the number of initial and retained cells after filtering, the reasons for cell removal, and a breakdown of the sample distribution.

The report includes:
- **Initial and retained cells**: The total number of cells before and after QC, along with the percentage of cells retained.
- **QC thresholds**: The thresholds used for gene count, UMI count, mitochondrial content, and hemoglobin content.
- **Filtered cells**: A table showing the number of cells removed for each QC reason.
- **Sample distribution**: A summary of how many cells remain in each sample after QC.

## Files Generated

1. **QC Plots**: Plots showing the distribution of QC metrics before and after filtering, saved as `QC_metrics_before_filtering.pdf` and `QC_metrics_after_filtering.pdf`.
2. **QC Report**: A text file (`QC_report.txt`) summarizing the QC steps and results.
3. **Processed Data**: The final processed Seurat object is saved as `scRNA_IPF_QC_Harmony_processed.RDS` or `scRNA_LCA_QC_Harmony_processed.RDS`.

## Usage

To run this pipeline, ensure that you have the required libraries installed and that the raw data files are properly placed in the respective directories.

1. **Install required R packages**:
   ```R
   install.packages(c("Seurat", "tidyverse", "patchwork", "harmony", "ggplot2"))
