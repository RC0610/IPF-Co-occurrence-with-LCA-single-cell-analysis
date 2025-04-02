####### EMT gene Analysis ########

# Step 1: Download EMT-related genes from GeneCards (https://www.genecards.org/). 
# Ensure to save the downloaded file as "EMT.csv". this data was downloaded at 2024.10.15.

#### 设置工作目录，建议使用相对路径
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Trajectory_and_EMT_Analysis")  # 请根据实际情况调整相对路径

#### 加载所需的R包
library(data.table)

#### 加载EMT.csv文件（确保文件路径正确）
EMT <- fread("EMT.csv")

#### 检查数据的前几行，确保加载正确
head(EMT)

#### 保存EMT相关基因为文本文件（相对路径）
write.table(EMT, file = "EMT_gene_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

