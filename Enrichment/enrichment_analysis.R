####################富集分析####################
####清空环境，设置报错语言
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
####2、加载R包
library(dplyr)
library(tidyverse)
library(DESeq2)
library(msigdbr)
library(dplyr)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(tidyverse)
library(data.table)
library(ggplot2)
setwd("H:\生信\IPF-Co-occurrence-with-LCA-single-cell-analysis\Enrichment")
####准备后续分析的基因名
# Step 1: Load the data files
IPF_time_diff_sig <- fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis/IPF_Time_diff_sig.csv")   # 在IPF中基底细胞状态的拟时间差异
IPF_BEAM_res_1 <- fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis/IPF_BEAM_res_1.csv")         # 在IPF中基底细胞的轨迹分支节点1差异
IPF_train_monocle_DEG <- read.table("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis/IPF.train.monocle.DEG.xls", header = TRUE, sep = "\t") # 在IPF中基底细胞状态相关的差异
EMT_genes <- fread("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis/EMT_gene_results.txt")             # EMT基因
filtered_DEGs <- read.table("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis/filtered.Degs.IPF.txt", header = TRUE, sep = "\t") # IPF的celltype差异基因
# Step 2: Extract gene columns based on their respective names
# Extract gene names from each file
genes_IPF_time_diff_sig <- IPF_time_diff_sig$gene_short_name
genes_IPF_BEAM_res_1 <- IPF_BEAM_res_1$gene_short_name
genes_IPF_train_monocle <- IPF_train_monocle_DEG$gene_short_name
genes_EMT <- EMT_genes$Gene Symbol
genes_filtered_DEGs <- filtered_DEGs$gene
genes_filtered_DEGs <- filtered_DEGs[filtered_DEGs$cluster == "Basal cell_IPF", "gene"]  # Modify column name if needed
# Step 3: Find the intersection of all gene lists
common_genes <- Reduce(intersect, list(genes_IPF_time_diff_sig, genes_IPF_BEAM_res_1, genes_IPF_train_monocle, genes_EMT, genes_filtered_DEGs))
####GO BP富集分析
#加载R包
library(org.Hs.eg.db)#人是Hs，鼠是Mm，需要安装
library(clusterProfiler)
erich.go.BP = enrichGO(gene =degs.list,#gene symbol
                       OrgDb = org.Hs.eg.db,#将symbol通过库转换
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,#p大于0.05删除
                       qvalueCutoff = 0.05)#同上，此处是GO通路的富集分析，结果erich.go.BP是s4对象，结果在里面的result里面，是一个data.frame的结构
#查看GO BP结果并输出保存为txt文档
erich.go.BP=erich.go.BP@result#查看通路富集的结果
write.table(erich.go.BP,"enrichment_results_GO_BP_Basal_IPF.txt",sep = "\t",col.names = NA)#文件名根据实际情况对应修改
#自定义作图GO BP
k = data.frame(erich.go.BP)
#作图的generatio是小数，所以要对generatio进行转换
x=k$BgRatio
x1=as.character(x)
a=strsplit(x1,split="/",fixed=T)#strsplit切割函数，主要切割字符
be=sapply(a,function(x){x[1]})#sapply函数提取before
be1=as.numeric(be)#转换为数值型
after=sapply(a,function(x){x[2]})#提取after
after=as.numeric(after)
k$GeneRatio = be1/after#geneRatio进行替换
#画图，用通道函数完成
font.size =13#字体大小
p=k %>% 
  ## 对进行p值排序
  arrange(pvalue) %>% 
  ##指定富集的通路数目
  slice(1:15) %>% 
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=pvalue, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

ggsave("GO_BP_Basal_IPF.pdf", p, width=12 ,height=8)  ##文件名根据实际情况对应修改
####kegg分析
#gene symbol转为ENTREZID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")#使用mapIDs把symbol转换为ENTREZID
#KEGG分析及作图气泡图
erich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "hsa",
                             keyType = "kegg")
#KEGG结果查看，并输出保存txt
zzh=erich.kegg.res@result#查看result结果
write.table(zzh,"enrichment_results_KEGG_Basal_IPF.txt",sep = "\t",col.names = NA)##文件名根据实际情况对应修改
#自定义作图KEGG
k = data.frame(zzh)#ggplot只能对data.frame进行作图
#作图的generatio是小数，所以要对generatio进行转换
x=k$GeneRatio
x1=as.character(x)
a=strsplit(x1,split="/",fixed=T)#strsplit切割函数，主要切割字符
be=sapply(a,function(x){x[1]})#sapply函数提取before
be1=as.numeric(be)#转换为数值型
after=sapply(a,function(x){x[2]})#提取after
after=as.numeric(after)
k$GeneRatio = be1/after#geneRatio进行替换
#画图，用通道函数完成
font.size =15#字体大小
p=k %>% 
  ## 对进行p值排序
  arrange(pvalue) %>% 
  ##指定富集的通路数目
  slice(1:7) %>% 
  ## 开始ggplot2 作图，其中fct_reorder调整因子level的顺序，10取决于对多少个通路作图
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  ## 画出点图
  geom_point(aes(color=pvalue, size = Count)) +
  ## 调整颜色，guide_colorbar调整色图的方向
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  ## 调整泡泡的大小
  scale_size_continuous(range=c(3, 8))+
  ## 如果用ylab("")或出现左侧空白
  labs(y=NULL) +
  ## 如果没有这一句，上方会到顶
  ggtitle("")+
  ## 设定主题
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))

ggsave("KEGG_Basal_IPF.pdf", p, width=12 ,height=8) #文件名根据实际情况对应修改
