####################IPF差异基因分析####################
####1、清空环境，设置报错语言
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
####2、加载所需R包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(SingleR)
library(data.table)
library(cowplot)
library(ggplot2)
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Differential_Expression_Analysis")
####3、前期质控、cluster的数据准备
scRNA=readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_IPF_harmony_clustered.RDS")##若为LCA，则加载scRNA_LCA_harmony_clustered.RDS
####4、根据文献标记基因手动注释 cluster################
ann=fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation/ann_IPF.csv") ##若为LCA，则改为ann_CA.csv
celltype = as.data.frame(ann)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
Idents(scRNA)="celltype"
#差异分析
Degs <- FindAllMarkers(object = scRNA, test.use = "wilcox", 
                          only.pos = TRUE, 
                          logfc.threshold = 0.25, 
                          p.adjust.method = "fdr")  # 使用 Benjamini-Hochberg 方法校正

# 筛选padj小于0.05、logFC大于0.25的marker基因
filtered.Degs <- Degs %>% 
  dplyr::select(gene, cluster, p_val, avg_log2FC, p_val_adj, pct.1, pct.2) %>% 
  subset(p_val_adj < 0.05 & avg_log2FC > 0.25 & pct.1 > 0.25 & pct.2 < 0.1)

# 保存筛选后的差异基因
write.table(filtered.Degs, file = "filtered.Degs.IPF.txt", sep = "\t", row.names = FALSE, quote = FALSE)##若为LCA时更换文件名为filtered.Degs.LCA.txt
