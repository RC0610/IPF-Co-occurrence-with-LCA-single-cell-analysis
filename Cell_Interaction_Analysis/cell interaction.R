####清空环境，设置报错语言
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
####加载R包
library(Seurat)
library(dplyr)
library(tidyverse)
library(data.table)
library(CellChat)
library(patchwork)
library(ggalluvial)
library(ggsci)
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Cell_Interaction_Analysis")
##读取工作文件
scRNA=readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_IPF_harmony_clustered.RDS")#若为LCA，则更换为scRNA_LCA_harmony_clustered.RDS
load("ref_Human_all.RData")#加载 CellChat 提供的人类 PPI 参考数据库
##细胞注释
celltype=fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation/ann_IPF.csv")#若为LCA，则更换为ann_CA.csv
celltype = as.data.frame(celltype)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
Idents(scRNA)="celltype"
##创建cellchat对象,提取表达矩阵和细胞分类信息
data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data")
identity <- subset(scRNA@meta.data, select = "celltype")#select选择针对什么进行细胞互作
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")#创建cellchat对象
##加载 CellChat 的人类数据库
CellChatDB <- CellChatDB.human
##对表达数据进行预处理
#将信号基因的表达数据进行子集化，以节省计算成本
cellchat <- subsetData(cellchat)
# 识别过表达基因,过滤表达低的基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
##相互作用推断
## 1、计算通信概率推断细胞互作的通信网络，得到相互作用的概率
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
## 2、如果特定细胞群中只有少数细胞，则过滤掉细胞间的通信，筛选出重要的互作，min.cells=3提示至少要在3个细胞中有互作才能留下
cellchat <- filterCommunication(cellchat, min.cells = 3)
#提取数据框，目的是找到自己感兴趣的通路
##在信号通路水平上推断细胞间的通讯，即细胞间的通讯网络
cellchat <- computeCommunProbPathway(cellchat)
##汇总通信概率来计算细胞间的聚合通信网络    
cellchat <- aggregateNet(cellchat)
## 3、计算聚合细胞互作通信网络
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
weight_matrix <- cellchat@net$weight
write.table(weight_matrix, file = "interaction_results_IPF.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)#根据具体情况修改为对应的文件名
#作图
# 加载R包
library(circlize)
library(dplyr)
#准备数据
links <- as.data.frame(as.table(weight_matrix))
colnames(links) <- c("from", "to", "weight")
links <- links %>% filter(weight > 0.05)
#设置 PDF 输出文件
library(RColorBrewer)
pdf("interaction_plots_IPF.pdf", width = 15, height = 15)  # 根据具体情况设置对应文件名
# Clear previous settings
circos.clear()
circos.par(start.degree = 90, gap.degree = 5)
# Define a color palette (choose a set of colors that fit your data)
colors <- brewer.pal(n = length(unique(c(links$from, links$to))), name = "Paired")
# Draw the chord diagram with the new color palette
chordDiagram(
  links,
  transparency = 0.5,
  grid.col = colors,
  annotationTrack = "grid",
  preAllocateTracks = 1
)
# Add cell type labels
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  },
  bg.border = NA
)
# Add title
title(main = "IPF cellchat", 
      cex.main = 1.5,
      font.main = 2,
      col.main = "black")
# 根据具体情况设置title
dev.off()
# 在cellchat中筛选wnt通路
wnt_pathway <- subsetCommunication(cellchat, pathways = "WNT")
# 转换通信数据为绘图格式
wnt_links <- as.data.frame(as.table(wnt_pathway$weight)) 
colnames(wnt_links) <- c("from", "to", "weight")
#过滤信号
wnt_links <- wnt_links %>% filter(weight > 0.05)
# 作图
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
# 获取Set2 调色板
set2_colors <- brewer.pal(8, "Paired")
# 可视化，并保持
pdf("wnt_in_IPF.pdf", width = 8, height = 8)#根据具体情况设置对应文件名
ggplot(data = wnt_links, aes(axis1 = from, axis2 = to, y = weight)) +
  geom_alluvium(aes(fill = from), width = 1/12) +  # Set width of alluvium
  geom_stratum(aes(fill = after_stat(stratum)), color = "black") +  # Set stratum color
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +  # Label the strata
  scale_x_discrete(limits = c("Source", "Target"), expand = c(0.05, 0.05)) +  # Label axes
  scale_fill_manual(values = set2_colors) +  # Apply color palette
  labs(title = "WNT in IPF",
       x = "Cell Types",
       y = "Communication Strength") +
  theme_minimal() +
  theme(legend.position = "none")  # Hide legend

dev.off()  # 根据具体情况设置对应title

