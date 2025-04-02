#### 1. 清空环境，设置报错语言
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

#### 2. 加载R包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(data.table)

setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation")

# 读取降维聚类后的数据
scRNA = readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_IPF_harmony_clustered.RDS")###做LCA时，将文件名替换为scRNA_LCA_harmony_clustered.RDS

########## 查找cluster marker gene并进行细胞注释 ############

## 使用 FindAllMarkers 函数查找每个 cluster 的marker基因
markers <- FindAllMarkers(object = scRNA, test.use = "wilcox", 
                          only.pos = TRUE, 
                          logfc.threshold = 0.25, 
                          p.adjust.method = "fdr")  # 使用 Benjamini-Hochberg 方法校正

# 筛选padj小于0.05、logFC大于0.25的marker基因
filtered.markers <- markers %>% 
  dplyr::select(gene, cluster, p_val, avg_log2FC, p_val_adj, pct.1, pct.2) %>% 
  subset(p_val_adj < 0.05 & avg_log2FC > 0.25 & pct.1 > 0.25 & pct.2 < 0.1)

# 保存筛选后的marker基因
write.table(filtered.markers, file = "filtered_IPF_marker_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)##若为LCA时更换文件名为filtered_LCA_marker_genes.txt
###The identified marker genes were then used to annotate cell clusters, and the cell groups were further validated by comparison to well-established cellular markers obtained from the literature.

########根据文献标记基因手动注释 cluster################
ann=fread("ann_IPF.csv") ##若为LCA，则改为ann_LCA.csv
celltype = as.data.frame(ann)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
Idents(scRNA)="celltype"
scRNA1=scRNA
#自定义作图IPF的celltype
pt_size = 1 #定义点的大小
df <- scRNA1@reductions$umap@cell.embeddings#提取坐标信息，每个细胞的横坐标和纵坐标的位置
df <- cbind(df, as.data.frame(scRNA1@active.ident))#df和active.ident都来自scRNA1，所以可不用merge合并，直接用cbine合并，横坐标是一一对应的
colnames(df) <- c("x", "y", "ident")#重命名列名为x，y和ident，ident是聚类数
meta=scRNA1@meta.data
meta=meta[,c(1,8)]
df.m=merge(df,meta,by=0)
df.m <- df.m %>%
  group_by(celltype) %>%
  summarise(
    x = median(x),
    y = median(y)
  )#每种细胞类型的中位数，x和y轴的中位数，按celltypes分类的

mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#6E568C","#7A142C","#E0367A")
scRNA1$celltype=factor(scRNA1$celltype,
                       levels = c("Basal cell_IPF","Type II pneumocyte_IPF","Ciliated cell_IPF","secretory cell_IPF","SLC16A7+ cell_IPF","Ionocyte cell_IPF","Macrophage_IPF","Myofibroblast_IPF"))


p1=DimPlot(scRNA1, reduction = "umap", 
           cols = mycol[1:10],
           group.by ="celltype",label = T)
ggsave("annotation_plot_IPF.pdf", p1, width=12 ,height=8)
#自定义作图LCA的celltype
pt_size = 1 #定义点的大小
df <- scRNA1@reductions$umap@cell.embeddings#提取坐标信息，每个细胞的横坐标和纵坐标的位置
df <- cbind(df, as.data.frame(scRNA1@active.ident))#df和active.ident都来自scRNA1，所以可不用merge合并，直接用cbine合并，横坐标是一一对应的
colnames(df) <- c("x", "y", "ident")#重命名列名为x，y和ident，ident是聚类数
meta=scRNA1@meta.data
meta=meta[,c(1,8)]
df.m=merge(df,meta,by=0)
df.m <- df.m %>%
  group_by(celltype) %>%
  summarise(
    x = median(x),
    y = median(y)
  )#每种细胞类型的中位数，x和y轴的中位数，按celltypes分类的

mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D")
scRNA1$celltype=factor(scRNA1$celltype,
                       levels = c("Basal cell_CA","Type II pneumocyte_CA","Ciliated cell_CA","Secretory cell_CA","SLC16A7+ cell_CA","Brush cell_CA","Macrophage_CA","Myeloid dendritic cell_CA","B cell_CA","T cell_CA"))


p1=DimPlot(scRNA1, reduction = "umap", 
           cols = mycol[1:10],
           group.by ="celltype",label = T)
ggsave("annotation_plot_LCA.pdf", p1, width=12 ,height=8)
