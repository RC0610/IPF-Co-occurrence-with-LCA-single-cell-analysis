####### Trajectory Analysis ########

# Step 1: Set working directory with absolute path (adjust as necessary)
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis")  # Ensure this path is correct for your system

#### 加载所需的R包
library(Seurat)          # 用于单细胞RNA数据分析
library(tidyverse)       # 包括dplyr等数据处理工具
library(dplyr)           # 数据清理和变换
library(patchwork)       # 用于组合ggplot2图形
library(harmony)         # 用于批次效应修正
library(data.table)      # 加载和处理大数据集
library(monocle)         # 用于轨迹分析

#### 数据准备，提取基底细胞
scRNA = readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_IPF_harmony_clustered.RDS")  # 更换为LCA数据时，路径应更改为scRNA_LCA_harmony_clustered.RDS
celltype = fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation/ann_IPF.csv")  # 更换为LCA时，路径应更改为ann_CA.csv

# Convert to data frame for processing
celltype = as.data.frame(celltype)

# Initialize celltype metadata as "NA"
scRNA@meta.data$celltype = "NA"

# Assign cell types based on cluster IDs in metadata
for(i in 1:nrow(celltype)) {
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]
}

# Subset the data to extract basal cells (change cluster name if different)
sc.t = scRNA[, rownames(subset(scRNA@meta.data, celltype == "Basal cell_IPF"))]  # Modify "Basal cell_IPF" based on your dataset

#### 创建 CellDataSet 对象
data <- as(as.matrix(sc.t@assays$RNA@counts), 'sparseMatrix')  # Create sparse matrix for monocle
pd <- new('AnnotatedDataFrame', data = sc.t@meta.data)  # Assign metadata
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))  # Assign gene information
fd <- new('AnnotatedDataFrame', data = fData)  # Assign feature data

# Create the CellDataSet (CDS) object for monocle
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())  # Create CDS object

#### 估计 size factor 和离散度
monocle_cds <- estimateSizeFactors(monocle_cds)  # Estimate normalization factors
monocle_cds <- estimateDispersions(monocle_cds)  # Estimate dispersions (variance in gene expression)

#### 过滤低质量的细胞
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)  # Detect genes with a minimum expression threshold

# 选择用于轨迹分析的基因（高变基因）
HSMM = monocle_cds  # Assign object to HSMM for convenience
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id  # Filter genes with high variance
HSMM <- setOrderingFilter(HSMM, disp.genes)  # Set ordering filter for trajectory
plot_ordering_genes(HSMM)  # Visualize the ordered genes

# Check the genes used for ordering
table(HSMM@featureData@data[["use_for_ordering"]])

#### 降维（选择方法DDRTree）
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')  # Reduce dimensions using DDRTree method, set to 2 components

# 根据轨迹顺序排列细胞并可视化
HSMM <- orderCells(HSMM)  # Order cells by trajectory

# Create the trajectory plots
p1 = plot_cell_trajectory(HSMM, color_by = "Pseudotime", size = 1, show_backbone = TRUE)  # Color by pseudotime
p2 = plot_cell_trajectory(HSMM, color_by = "State", size = 1, show_backbone = TRUE)  # Color by cell state/celltype
p = p1 + p2  # Combine the plots

p4 = plot_cell_trajectory(HSMM, color_by = "Trajectory") + facet_wrap("~Trajectory", nrow = 1)  # Split trajectory plot by cell type
ggsave("IPF_Pseudotime_celltype.pdf", plot = p1 | p2, width = 12, height = 6)  # Save the plot
ggsave("IPF_trajectory_facet.pdf", plot = p4, width = 12, height = 6)  # Save the facet plot

#### 轨迹差异分析
diff = differentialGeneTest(HSMM[disp.genes,], fullModelFormulaStr = "~State", cores = 1)  # 基于细胞状态的差异基因检测，若为celltype则为基因细胞类型的差异基因检测
deg = subset(diff, qval < 0.01)  # Filter significant genes (q-value < 0.01)
deg = deg[order(deg$qval, decreasing = F), ]  # Order by q-value
write.table(deg, file = "IPF.train.monocle.DEG.xls", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)  # Save the DEG results

# Extract the genes of interest
ordergene = rownames(deg)  # Get the gene names of the ordered genes
cds_subset = HSMM[c("KRT5", "PTTG1", "S100A14", "RND3", "SERPINB3", "KRT6A", "KRT17", "IGFBP5", "KRT19", "TAGLN"), ]  # Subset of genes related with EMT, Trajectory, and DEGs.

# Plot the genes in pseudotime
p1 = plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 = plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
p = p1 | p2  # Combine the plots

ggsave("gene.trajectory.pdf", plot = p, width = 12, height = 9)  # Save the plot

#### 拟时差异基因分析
Time_diff = differentialGeneTest(HSMM[ordergene,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")  # Differential analysis based on pseudotime
write.csv(Time_diff, "IPF_Time_diff_all.csv", row.names = F)  # Save all pseudotime-related genes

Time_genes = as.character(c("KRT5", "PTTG1", "S100A14", "RND3", "SERPINB3", "KRT6A", "KRT17", "IGFBP5", "KRT19", "TAGLN"))  # Subset of genes related with EMT, Trajectory, and DEGs.
p = plot_pseudotime_heatmap(HSMM[Time_genes,], num_clusters = 4, show_rownames = TRUE, return_heatmap = TRUE)  # Plot heatmap
ggsave("IPF_plot_pseudotime_heatmap.pdf", plot = p, width = 12, height = 6)  # Save the heatmap

# Sort the genes by heatmap results and save the file
hp.genes = p$tree_row$labels[p$tree_row$order]
Time_diff_sig = Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "IPF_Time_diff_sig.csv", row.names = F)  # Save significant pseudotime genes

#### 分支分析（BEAM）
marker_genes = row.names(subset(fData(HSMM), gene_short_name %in% c("KRT5", "PTTG1", "S100A14", "RND3", "SERPINB3", "KRT6A", "KRT17", "IGFBP5", "KRT19", "TAGLN")))  # Define marker genes for BEAM
diff_test_res = differentialGeneTest(HSMM[marker_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)")  # Perform differential gene test for BEAM

#### 单细胞轨迹的分支分析（Branch analysis）
plot_cell_trajectory(HSMM, color_by = "State")  # Plot the cell trajectory colored by State

# BEAM statistical analysis
BEAM_res = BEAM(HSMM[ordergene,], branch_point = 1, cores = 2, progenitor_method = 'duplicate')  # Perform BEAM for branch analysis，choosing branch_point by nodes of gene.trajectory.pdf
BEAM_res = BEAM_res[order(BEAM_res$qval), ]  # Order by q-value
BEAM_res = BEAM_res[, c("gene_short_name", "pval", "qval")]  # Select relevant columns
BEAM_res = BEAM_res[row.names(subset(BEAM_res, qval < 0.01)), ]  # Filter based on q-value
head(BEAM_res)

# Save BEAM results
write.csv(BEAM_res, "IPF_BEAM_res_1.csv", row.names = F)  # Save the BEAM results,file name was changed by branch_point

# Plot BEAM heatmap for selected genes, file name was changed by branch_point
pdf("IPF_BEAM_heatmap_1.pdf", width = 12, height = 10)
plot_genes_branched_heatmap(HSMM[Time_genes,], branch_point = 1, num_clusters = 4, cores = 2, use_gene_short_name = TRUE, show_rownames = TRUE)
dev.off()  # Close the PDF
