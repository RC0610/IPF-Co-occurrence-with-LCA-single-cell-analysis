###########SCENIC##############
####清空环境，设置报错语言
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
####加载R包
library(sp)
library(SeuratObject)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(harmony)
library(data.table)
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/SCENIC_Analysis")
##加载源分析文件
scRNA=readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Clustering_Analysis/scRNA_IPF_harmony_clustered.RDS")#若为LCA，则更换为scRNA_LCA_harmony_clustered.RDS
##细胞注释
celltype=fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Marker_gene_Expression_Analysis_and_Cell_Annotation/ann_IPF.csv")#若为LCA，则更换为ann_CA.csv
celltype= as.data.frame(celltype)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
scRNA=scRNA[,rownames(subset(scRNA@meta.data,celltype== "Basal cell_IPF"))]#根据所需要分析的细胞进行修改

subcell.id <- sample(colnames(scRNA),1500) #抽样1500个细胞进行scenic的分析
scRNAsub <- scRNA[,subcell.id]
exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
##设置分析环境
mydbDIR <- "H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/SCENIC_Analysis" #设置数据存放的路径
mydbs <- c( "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
names(mydbs) <- c("10kb","500bp")
##加载数据资源
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
#初始化 SCENIC 设置,设置分析环境
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=4,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "os")
saveRDS(scenicOptions, "int/scenicOptions.rds")
##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达，此步骤过滤后需要检查,检查关注的基因是否被过滤掉
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]

##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
#根据表达数据推断潜在的转录因子靶标，使用Genie3推断转录因子与靶基因之间的关系
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
###共表达网络直接载入之前的集群结果
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #1. 获取共表达模块，需要改为int，只能识别int下的结果
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)  #2. 获取regulons，这里的h38_10kb和h38_100kb需要在mydbDIR设置的这个工作环境之下
library(foreach)
exprMat_all <- as.matrix(scRNA@assays$RNA@counts)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all) 
aucell <- readRDS("int/3.4_regulonAUC.Rds")
auc_mtx=aucell@assays@data@listData$AUC
auc_mtx=t(auc_mtx)
# 获取 AUC 矩阵和基因表达矩阵
auc_df <- as.data.frame(auc_mtx)
# 提取基因表达矩阵（假设 Seurat 对象中基因表达数据在 `RNA` assay 中）
gene_expr <- as.data.frame(GetAssayData(scRNA, slot = "data"))
# 确保 AUC 和基因表达矩阵有相同的细胞顺序
common_cells <- intersect(rownames(auc_df), colnames(gene_expr))
grep("JUN",colnames(auc_df),value = T)#使用grep进行搜索，搜索含有JUN的通路名字

auc_df <- auc_df[common_cells,c("ATF3 (4182g)","BCLAF1 (4277g)","EBF1 (28g)","EGR1 (4742g)","NFATC4 (61g)","JUND (240g)","KLF5 (441g)","TFAP2A (110g)","SOX9_extended (18g)","EHF (2417g)","ATF4 (103g)","EGR3 (122g)","FOSL1 (240g)","JUN (95g)") ]
gene=gene_expr[c("KRT5","PTTG1","S100A14","RND3","SERPINB3","KRT6A","KRT17","IGFBP5","KRT19","TAGLN"),common_cells]
# 计算每个基因与每个 TF AUC 之间的相关性
cor_matrix <- cor(t(gene), auc_df, method = "pearson")

# 绘制相关性热图
pdf("scenic_plots_IPF_Basal.pdf",width = 12,height = 8)
pheatmap(cor_matrix,
         scale = "none",
         cluster_rows = TRUE, # 基因聚类
         cluster_cols = TRUE, # TF 聚类
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Gene Expression and TF AUC Correlation")

dev.off()#文件名根据实际情况对应修改
