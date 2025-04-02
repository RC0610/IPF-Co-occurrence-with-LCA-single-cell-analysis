#################### 降维与聚类 ####################
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis")
library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)

#### 1. 加载质控后数据 ####
scRNA <- readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Preprocessing_and_QC/scRNA_IPF_QC_Harmony_processed.RDS")#若为LCA，则改为scRNA_LCA_QC_Harmony_processed.RDS

#### 2. PCA降维 ####
# 使用 PCA + ElbowPlot 选择前 20 个主成分
# 分辨率 0.5 是 Seurat 的推荐初始值
scRNA <- RunPCA(scRNA, npcs = 50, verbose = TRUE)

# 可视化判断最佳维度数
pdf("ElbowPlot.pdf", width = 6, height = 4)
dev.off()

#### 3. Harmony批次校正（数据质控已经完成，故跳过）####

#### 4. 聚类分析 ####
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 0.5)  # 可根据分辨率优化

#### 5. 非线性降维可视化 ####
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20)

#### 6. 可视化结果 ####
p <- DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 1.5)

pdf("DimPlot_UMAP.pdf", width = 12, height = 6)
p
dev.off()
#### 7. 保存处理结果 ####
saveRDS(scRNA, file = "scRNA_IPF_harmony_clustered.RDS")##若为LCA，则改为scRNA_LCA_harmony_clustered.RDS
#### 8. 保存当前R环境信息####
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

clustering_results <- data.frame(
  Cluster = names(table(scRNA$seurat_clusters)),  # 聚类编号
  Cell_Count = as.vector(table(scRNA$seurat_clusters))  # 每个聚类的细胞数
)

# 保存聚类结果为文本文件
write.table(clustering_results, file = "clustering_results_IPF.txt", sep = "\t", row.names = FALSE, quote = FALSE)##若为LCA，则改为clustering_results_LCA.txt
