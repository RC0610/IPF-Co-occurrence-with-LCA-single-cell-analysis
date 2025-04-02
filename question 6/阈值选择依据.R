##加载r包
library(magrittr)  
library(tidyverse)  
library(Seurat)
library(patchwork)
setwd("C:/Users/Lenovo/Desktop/1")

#### 选择阈值
scRNA=readRDS("scRNA.CA.RDS") #质控前的scRNA，即读取的原始数据
quantile(scRNA$nCount_RNA, probs = c(0.95, 0.99, 0.995, 1.00))

###通过nfeature和ncount的线性关系证明高 nCount_RNA 的细胞是高质量、非异常的
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  #检查线性关，质控前的scRNA，即读取的原始数据

####检查这些细胞是否为特定亚群贡献主力（高表达细胞分布广泛，并不集中于某 cluster，那说明没有引起严重偏差）scRNA为经过质控后PCA、UMAP、聚类后的scRNA
scRNA=readRDS("scRNA_LCA_harmony_clustered.RDS")
high_ncount_cells <- WhichCells(scRNA, expression = nCount_RNA > 50000 & nCount_RNA <= 200000)
length(high_ncount_cells)  # 看数量
DimPlot(scRNA, group.by = "seurat_clusters", cells.highlight = high_ncount_cells, cols.highlight = "red") # 可视化这批细胞在 UMAP 上的分布

###对比不同阈值下的主要结论一致性
# 原始对象，读取的原始scRNA
scRNA_raw <- readRDS("scRNA.CA.RDS")
###质控指标计算 ####
scRNA_raw[["percent.mt"]] <- PercentageFeatureSet(scRNA_raw, pattern = "^MT-")
hb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
hb_genes <- intersect(hb_genes, rownames(scRNA_raw))
scRNA_raw[["percent.HB"]] <- PercentageFeatureSet(scRNA_raw, features = hb_genes)

# A. 使用 50,000 阈值过滤
scRNA_50K <- subset(scRNA_raw, subset =
                      nFeature_RNA > 300 & nFeature_RNA < 9000 &
                      nCount_RNA > 500 & nCount_RNA < 50000 &
                      percent.mt < 20 & percent.HB < 5
)

# B. 使用 200,000 阈值过滤（当前设置）
scRNA_200K <- subset(scRNA_raw, subset =
                       nFeature_RNA > 300 & nFeature_RNA < 9000 &
                       nCount_RNA > 500 & nCount_RNA < 200000 &
                       percent.mt < 20 & percent.HB < 5
)
###对两组数据分别进行标准化、PCA、UMAP、聚类
process_pipeline <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj) %>%
    FindVariableFeatures(nfeatures = 3000) %>%
    ScaleData(vars.to.regress = c("percent.mt")) %>%
    RunPCA(npcs = 50) %>%
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.5)
  return(seurat_obj)
}

scRNA_50K <- process_pipeline(scRNA_50K)
scRNA_200K <- process_pipeline(scRNA_200K)
###UMAP对比
p1 <- DimPlot(scRNA_50K, reduction = "umap", label = TRUE) + ggtitle("nCount < 50,000")
p2 <- DimPlot(scRNA_200K, reduction = "umap", label = TRUE) + ggtitle("nCount < 200,000")

pdf("Compare_UMAP_thresholds.pdf", width = 10, height = 5)
p1 + p2
dev.off()
###Cluster组成和数目对比
table_50K <- table(scRNA_50K$seurat_clusters)
table_200K <- table(scRNA_200K$seurat_clusters)

print("50K Clusters:")
print(table_50K)

print("200K Clusters:")
print(table_200K)
###关键 marker 鉴定与表达趋势比较
p1 <- FeaturePlot(scRNA_50K, features = c("KRT5"), reduction = "umap") + 
  ggtitle("KRT5 (Basal cell marker) - 50K threshold")

p2 <- FeaturePlot(scRNA_200K, features = c("KRT5"), reduction = "umap") + 
  ggtitle("KRT5 (Basal cell marker) - 200K threshold")

p3 <- FeaturePlot(scRNA_50K, features = c("SFTPC"), reduction = "umap") + 
  ggtitle("SFTPC (AT II marker) - 50K threshold")

p4 <- FeaturePlot(scRNA_200K, features = c("SFTPC"), reduction = "umap") + 
  ggtitle("SFTPC (AT II marker) - 200K threshold")

# 拼接图像：上下对比，左右排列
combined_plot <- (p1 | p2) / (p3 | p4)

# 展示图像
print(combined_plot)

# 可选择保存
ggsave("FeaturePlot_KRT5_SFTPC_compare_50K_vs_200K.pdf", combined_plot, width = 12, height = 10)
