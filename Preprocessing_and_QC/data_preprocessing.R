#################### 数据读入 ####################
#### 1. 初始化环境 ####
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

#### 2. 设置工作路径 ####
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Preprocessing_and_QC")
data_dir <- file.path("data")
h5_dir <- file.path(data_dir, "h5")       # GSE143706（需包含5个h5文件）
matrix1_dir <- file.path(data_dir, "matrix1")  # GSE128033（需包含3个子目录）
matrix2_dir <- file.path(data_dir, "matrix2")  # GSE190889（需包含5个子目录）

#### 3. 包加载 ####
library(Seurat)
library(tidyverse)
library(patchwork)

#### 4. 自定义函数（新增样本重命名逻辑）####
read_sc_h5 <- function(h5_dir, dataset_id, start_sample) {
  if (!dir.exists(h5_dir)) stop("h5目录不存在：", normalizePath(h5_dir))
  
  h5_files <- list.files(h5_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
  if (length(h5_files) == 0) stop("未找到.h5文件")
  
  message("读取h5文件：\n", paste(h5_files, collapse = "\n"))
  
  # 生成样本编号 sample[start_sample] 到 sample[start_sample + n -1]
  sample_ids <- sprintf("sample%d", seq(start_sample, start_sample + length(h5_files) - 1))
  
  sc_list <- lapply(seq_along(h5_files), function(i) {
    f <- h5_files[i]
    counts <- Read10X_h5(f)
    rownames(counts) <- gsub("_", "-", rownames(counts))
    
    obj <- CreateSeuratObject(
      counts,
      min.cells = 3,
      min.features = 300,
      project = dataset_id  # orig.ident设为数据集ID
    )
    obj$sample_id <- sample_ids[i]  # 使用生成的样本编号
    obj$dataset <- dataset_id
    obj
  })
  
  merge(sc_list[[1]], sc_list[-1])
}

read_sc_matrix <- function(matrix_dir, dataset_id, start_sample) {
  if (!dir.exists(matrix_dir)) stop("目录不存在：", normalizePath(matrix_dir))
  
  samples <- list.dirs(matrix_dir, full.names = FALSE, recursive = FALSE)
  if (length(samples) == 0) stop("未找到样本子目录")
  
  message("处理样本：\n", paste(samples, collapse = "\n"))
  
  # 生成样本编号 sample[start_sample] 到 sample[start_sample + n -1]
  sample_ids <- sprintf("sample%d", seq(start_sample, start_sample + length(samples) - 1))
  
  sc_list <- lapply(seq_along(samples), function(i) {
    s <- samples[i]
    data_path <- file.path(matrix_dir, s)
    counts <- Read10X(data_path)
    rownames(counts) <- gsub("_", "-", rownames(counts))
    
    obj <- CreateSeuratObject(
      counts,
      min.cells = 3,
      min.features = 300,
      project = dataset_id
    )
    obj$sample_id <- sample_ids[i]  # 使用生成的样本编号
    obj$dataset <- dataset_id
    obj
  })
  
  merge(sc_list[[1]], sc_list[-1])
}

#### 5. 数据读取主流程（指定样本起始编号）####
# GSE143706 (5个样本) → sample1-sample5
scRNA_h5 <- read_sc_h5(h5_dir, dataset_id = "GSE143706", start_sample = 1)

# GSE128033 (3个样本) → sample6-sample8
scRNA_mat1 <- read_sc_matrix(matrix1_dir, dataset_id = "GSE128033", start_sample = 6)

# GSE190889 (5个样本) → sample9-sample13
scRNA_mat2 <- read_sc_matrix(matrix2_dir, dataset_id = "GSE190889", start_sample = 11)

#### 6. 数据合并与验证 ####
scRNA <- merge(scRNA_h5, y = c(scRNA_mat1, scRNA_mat2))

# 验证样本命名
meta_data <- scRNA@meta.data
cat("样本分布验证:\n")
print(table(meta_data$dataset, meta_data$sample_id))

#### 7. 数据保存 ####
saveRDS(scRNA, file = "scRNA.IPF.RDS")##若为LCA数据集，则此处为scRNA.LCA.RDS


####################数据质控####################
#### 1. 初始化环境 ####
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Preprocessing_and_QC")
#### 2. 包加载 ####
library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)
library(knitr)
#### 3. 数据加载 ####
scRNA <- readRDS("H:/生信/IPF-Co-occurrence-with-LCA-single-cell-analysis/Preprocessing_and_QC/scRNA.IPF.RDS")## 如果是LCA数据，可以修改为："scRNA.LCA.RDS"
ncol_original <- ncol(scRNA)
#### 4. 统一质控标准 ####
QC_THRESHOLDS <- list(
  nFeature_RNA = c(300, 7500),     # 基因数/细胞: 
  # - 下限300: 排除空液滴或低质量细胞 (参考 Satija Lab建议)
  # - 上限7500: 排除多细胞或高RNA含量的异常细胞 (基于数据分布的99%分位数)
  
  nCount_RNA = c(500, 50000),      # UMI数/细胞:
  # - 下限500: 确保最低测序深度 (参考 Satija Lab建议)
  # - 上限50000: 用于排除潜在的多细胞或过量RNA含量细胞（双细胞、细胞团的异常细胞）(基于数据分布的99.5%分位数)
  
  percent_mt = 20,                 # 线粒体基因百分比:
  # - 阈值20%: 超过此值提示细胞损伤 (参考 Satija Lab建议)
  
  percent_hb = 5                   # 血红蛋白基因百分比:
  # - 血红蛋白比例上限5%：用于去除可能来源于红细胞的污染（常用于肺组织等含血丰富组织）（PMID：26216973）
)

#### 5. 质控指标计算 ####
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
hb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
hb_genes <- intersect(hb_genes, rownames(scRNA))
scRNA[["percent.HB"]] <- PercentageFeatureSet(scRNA, features = hb_genes)

#### 6. 样本级信息提取（关键修复）####
print(table(scRNA$orig.ident))  # 应显示 sample1-sample13

#### 7. 质控可视化函数 ####
qc_plots <- function(obj, title) {
  plot_list <- list(
    VlnPlot(obj, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0) +
      geom_hline(yintercept = QC_THRESHOLDS$nFeature_RNA, linetype = "dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    
    VlnPlot(obj, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0) +
      geom_hline(yintercept = QC_THRESHOLDS$nCount_RNA, linetype = "dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    
    VlnPlot(obj, features = "percent.mt", group.by = "orig.ident", pt.size = 0) +
      geom_hline(yintercept = QC_THRESHOLDS$percent_mt, linetype = "dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    
    VlnPlot(obj, features = "percent.HB", group.by = "orig.ident", pt.size = 0) +
      geom_hline(yintercept = QC_THRESHOLDS$percent_hb, linetype = "dashed", color = "red") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  
  wrap_plots(plot_list, ncol = 2) + 
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
}

#### 8. 生成质控图表 ####
pdf("QC_metrics_before_filtering.pdf", width = 16, height = 12)
qc_plots(scRNA, "Pre-QC Metrics (Sample Level)")
dev.off()

#### 9. 数据过滤 ####
scRNA_pre <- scRNA
scRNA <- subset(scRNA, subset = 
                  nFeature_RNA > QC_THRESHOLDS$nFeature_RNA[1] & 
                  nFeature_RNA < QC_THRESHOLDS$nFeature_RNA[2] &
                  nCount_RNA > QC_THRESHOLDS$nCount_RNA[1] &
                  nCount_RNA < QC_THRESHOLDS$nCount_RNA[2] &
                  percent.mt < QC_THRESHOLDS$percent_mt &
                  percent.HB < QC_THRESHOLDS$percent_hb
)

pdf("QC_metrics_after_filtering.pdf", width = 16, height = 12)
qc_plots(scRNA, "Post-QC Metrics (Sample Level)")
dev.off()

#### 10. 细胞周期处理 ####
scRNA <- CellCycleScoring(scRNA,
                          s.features = intersect(cc.genes$s.genes, rownames(scRNA)),
                          g2m.features = intersect(cc.genes$g2m.genes, rownames(scRNA)))

#### 11. 标准化与批次校正 ####
scRNA <- NormalizeData(scRNA) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) %>%
  RunPCA(npcs = 50)

scRNA <- RunHarmony(scRNA, 
                    group.by.vars = "orig.ident",
                    theta = 2,
                    plot_convergence = TRUE)

#### 12. 质控验证图表（修复分组）####
p_int <- DimPlot(scRNA, reduction = "harmony", group.by = "orig.ident") +
  ggtitle("Batch Integration by Dataset")

p_smp <- DimPlot(scRNA, reduction = "harmony", group.by = "orig.ident") +  # 使用 orig.ident
  ggtitle("Batch Integration by Sample")

pdf("QC_integration_validation.pdf", width = 14, height = 6)
p_int + p_smp
dev.off()

#### 13. 生成质控报告 ####
batch_stats <- as.data.frame(table(Sample_ID = scRNA$orig.ident))
colnames(batch_stats) <- c("Sample_ID", "Cells_After_QC")

# 计算各阈值过滤的细胞数（来自 scRNA_pre）
filtered_stats <- data.frame(
  Reason = c("Low nFeature_RNA", "High nFeature_RNA", 
             "Low nCount_RNA", "High nCount_RNA",
             "High MT%", "High HB%"),
  Cells = c(
    sum(scRNA_pre$nFeature_RNA < QC_THRESHOLDS$nFeature_RNA[1]),
    sum(scRNA_pre$nFeature_RNA > QC_THRESHOLDS$nFeature_RNA[2]),
    sum(scRNA_pre$nCount_RNA < QC_THRESHOLDS$nCount_RNA[1]),
    sum(scRNA_pre$nCount_RNA > QC_THRESHOLDS$nCount_RNA[2]),
    sum(scRNA_pre$percent.mt >= QC_THRESHOLDS$percent_mt),
    sum(scRNA_pre$percent.HB >= QC_THRESHOLDS$percent_hb)
  )
)

qc_report <- paste(
  "=== QUALITY CONTROL REPORT ===",
  paste0("\nGenerated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("\nInitial cells: ", ncol_original),
  paste0("\nRetained cells: ", ncol(scRNA), " (", 
         round(ncol(scRNA)/ncol_original*100, 1), "%)"),
  
  "\n\n=== QC THRESHOLDS APPLIED ===",
  capture.output(kable(data.frame(
    Parameter = c("Genes per Cell (nFeature_RNA)", 
                  "UMIs per Cell (nCount_RNA)", 
                  "Mitochondrial Content (%)", 
                  "Hemoglobin Content (%)"),
    Min = c(QC_THRESHOLDS$nFeature_RNA[1], 
            QC_THRESHOLDS$nCount_RNA[1],
            NA, NA),
    Max = c(QC_THRESHOLDS$nFeature_RNA[2],
            QC_THRESHOLDS$nCount_RNA[2],
            QC_THRESHOLDS$percent_mt,
            QC_THRESHOLDS$percent_hb),
    Rationale = c(
      "300–7500: Filters empty droplets or stressed cells (Satija Lab; top 99%)",
      "500–50000: Ensures minimum sequencing depth; removes multiplets (>99.5%)",
      "<20%: High mitochondrial % suggests damaged/apoptotic cells",
      "<5%: High HB% may indicate erythrocyte contamination (PMID:26216973)"
    )
  ), align = 'c', caption = "Applied QC Thresholds")),
  
  "\n\n=== FILTERED CELLS BY CRITERIA ===",
  capture.output(kable(filtered_stats, align = 'c', 
                       col.names = c("Filter Reason", "Number of Cells Removed"))),
  
  "\n\n=== SAMPLE DISTRIBUTION ===",
  capture.output(kable(batch_stats, digits = 0,
                       col.names = c("Sample ID", "Cells After QC"))),
  
  "\n\n=== ADDITIONAL NOTES ===",
  "1. Thresholds validated via QC distribution plots (see PDF)",
  "2. Batch correction performed via Harmony (theta=2, 50 PCs)",
  "3. Cell cycle scores (S, G2M) regressed in scaling step",
  "4. QC plots saved as QC_metrics_[before/after]_filtering.pdf",
  sep = "\n"
)

writeLines(qc_report, "QC_report.txt")
#### 14. 保存批次校正后的数据 ####
saveRDS(scRNA, file = "scRNA_IPF_QC_Harmony_processed.RDS") ###若为LCA数据集，则此处为scRNA_LCA_QC_Harmony_processed.RDS
