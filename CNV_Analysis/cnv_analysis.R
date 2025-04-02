####清空环境，设置报错语言
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
####加载R包
library(dplyr)
library(Seurat)
library(tidyverse)
library(infercnv)
library(patchwork)
library(data.table)
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/CNV_Analysis")
##读取工作文件
scRNA=readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_LCA_harmony_clustered.RDS")
celltype=fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation/ann_CA.csv")
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
############################CNV分析###############################################
#BiocManager::install("infercnv")
library(infercnv)
pos=read.table("human.gene.positions")#读取包含基因位置信息的文件
pos1=distinct(pos,V7,.keep_all = TRUE)
rownames(pos1)=pos1$V7
pos2=select(pos1,V7,V2,V3,V4)
write.table(pos2, 'geneLocate.txt', row.names=F, col.names=F, sep='\t')
scRNA1=scRNA[,sample(colnames(scRNA),500)]
exprMatrix <- GetAssayData(scRNA1, slot = 'counts')
cellAnnota <- subset(scRNA1@meta.data, select='celltype')
groupFiles='groupFiles.txt'
dim(exprMatrix)
write.table(cellAnnota,file ="groupFiles.txt",sep = '\t',col.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprMatrix,
                                    annotations_file="groupFiles.txt",
                                    delim="\t",
                                    gene_order_file="geneLocate.txt",
                                    ref_group_names=c("B cell_CA","T cell_CA","Myeloid dendritic cell_CA"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir=  'cnv1/' ,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   #  cluster_by_groups：先区分细胞来源，再做层次聚类
                             k_obs_groups = 1,
                             HMM = T,
                             hclust_method="ward.D2", plot_steps=F,denoise=TRUE,
                             num_threads=4,
                             write_expr_matrix=T)
##可视化分组CNV打分
grp=read.table("cnv1/infercnv.observation_groupings.txt",sep = "",header = T)
obs=read.table("cnv1/infercnv.observations.txt", header = T,check.names = F)
max(obs)
min(obs)
obs[obs>0.9 & obs<0.93]=2
obs[obs>=0.93 & obs<0.95]=1
obs[obs>=0.95 & obs<1.05]=0
obs[obs>=1.05 & obs<1.07]=1
obs[obs>=1.07 & obs<1.1]=2

scores=as.data.frame(colSums(obs))
scores$celltype=grp$Dendrogram.Group
colnames(scores)=c("score","celltype_s")
scores=data.frame(scores)
scores$celltype <- gsub("_s.*", "", scores$celltype)
scores=scores[,c(1,3)]
scores$celltype=factor(scores$celltype,
                       levels = c("Type II pneumocyte_CA","Basal cell_CA","Macrophage_CA","Brush cell_CA","Ciliated cell_CA","Secretory cell_CA","SLC16A7+ cell_CA"))
##CNV打分及可视化
library(ggpubr)
p1=ggboxplot(scores,"celltype","score",fill = "celltype")+
  scale_fill_manual(values = c("Type II pneumocyte_CA" = "#D20A13", "Basal cell_CA" = "#223D6C", "Macrophage_CA" = "#7A142C","Brush cell_CA" = "#58CDD9","Ciliated cell_CA" = "#FFD121","Secretory cell_CA" = "#088247","SLC16A7+ cell_CA" ="#11AA4D")) +
  labs(fill = "cluster")
ggsave("cnv boxplot.pdf", plot = p1, width = 12, height = 6) 

##CNV的热图
library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "infercnv",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("blue","white","red"), c(2, 2))) #改颜色
