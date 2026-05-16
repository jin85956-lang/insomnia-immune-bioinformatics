# ============================================================
# 00_setup_and_data.R
# ============================================================
# 对应 Methods §2.1 Data acquisition and preprocessing
#
# 输入：GEO 数据库在线下载 GSE208668
# 输出：
#   Data/GSE208668_processed.RData  (expr_final, group, pheno)
#   Data/GSE208668_expr_matrix.csv  (14539 基因 × 42 样本)
#   Data/GSE208668_pheno.csv        (42 样本 × 6 列)
# ============================================================

# ---- 0. 项目根目录与子文件夹 ----
project_dir <- "F:/XiaoChen_Course/Insomnia_Project"
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

folders <- c(
  "Phase0_eQTL_check",
  "Phase1_GSVA",
  "Phase2_WGCNA",
  "Phase3_ML",
  "Phase3.5_CrossCohort",
  "Phase4_Immune",
  "Phase5_SMR",
  "Phase6_DGIdb",
  "Phase6_Writing",
  "Data",
  "Results"
)
for (f in folders) {
  dir.create(file.path(project_dir, f), showWarnings = FALSE)
}

setwd(project_dir)
cat("项目目录就绪：", project_dir, "\n")

# ---- 1. 加载依赖 ----
library(GEOquery)
library(limma)
library(Biobase)

# ---- 2. 下载 GSE208668 表达数据 ----
# GSE208668：Chronic insomnia PBMC microarray, Illumina HumanHT-12 V4 (GPL10904)
# 42 samples: 17 insomnia + 25 control (60-80 岁)
# 参考：Cho 2022, PMID 待补
cat("下载 GSE208668 ...\n")
gse <- getGEO("GSE208668", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]

# ---- 3. 提取原始表达矩阵 ----
expr_raw <- exprs(eset)
cat("原始探针数:", nrow(expr_raw), "\n")
cat("样本数:    ", ncol(expr_raw), "\n")

# ---- 4. 过滤无效探针 ----
# 4a. 删掉以 "HS." 或 "LOC" 开头的低质量探针
all_ids <- rownames(expr_raw)
keep_clean <- !grepl("^HS\\.|^LOC", all_ids)
expr_clean <- expr_raw[keep_clean, ]
# 4b. 删掉含 NA 的行
expr_clean <- expr_clean[complete.cases(expr_clean), ]
cat("清洗后基因数:", nrow(expr_clean), "\n")

# 4c. 去基因名重复（保留平均表达最高的探针）
if (any(duplicated(rownames(expr_clean)))) {
  avg <- rowMeans(expr_clean, na.rm = TRUE)
  o <- order(avg, decreasing = TRUE)
  expr_clean <- expr_clean[o, ]
  expr_clean <- expr_clean[!duplicated(rownames(expr_clean)), ]
}
cat("去重后基因数:", nrow(expr_clean), "\n")

# ---- 5. 构建分组变量 ----
# pData(eset)$`insomnia:ch1` 的取值: "yes"=Insomnia, "no"=Control
group <- factor(pData(eset)$`insomnia:ch1`, levels = c("no", "yes"))
cat("分组分布:\n")
print(table(group))

# ---- 6. 低表达过滤 ----
# 标准：基因在 insomnia 或 control 组中至少 50% 样本表达值 > 7.0（log2 scale）
insomnia_idx <- which(group == "yes")
control_idx  <- which(group == "no")

keep_expressed <- apply(expr_clean, 1, function(x) {
  (sum(x[insomnia_idx] > 7.0) >= length(insomnia_idx) * 0.5) |
  (sum(x[control_idx]  > 7.0) >= length(control_idx)  * 0.5)
})

expr_final <- expr_clean[keep_expressed, ]
cat("\n=== 最终表达矩阵 ===\n")
cat("基因数:", nrow(expr_final), "\n")   # 14,539
cat("样本数:", ncol(expr_final), "\n")   # 42

# ---- 7. 构建表型数据 ----
pheno <- data.frame(
  SampleID = colnames(expr_final),
  group    = as.character(group),
  age      = as.numeric(pData(eset)$`age:ch1`),
  gender   = pData(eset)$`gender:ch1`,
  bmi      = as.numeric(pData(eset)$`bmi:ch1`),
  bdi      = as.numeric(pData(eset)$`bdi:ch1`),
  stringsAsFactors = FALSE
)
cat("\n表型数据前 6 行:\n")
print(head(pheno))

# ---- 8. 保存 ----
save(expr_final, group, pheno,
     file = "Data/GSE208668_processed.RData")
write.csv(expr_final, "Data/GSE208668_expr_matrix.csv")
write.csv(pheno, "Data/GSE208668_pheno.csv", row.names = FALSE)

cat("\n========== 00 数据预处理完成 ==========\n")
cat("expr_final:", nrow(expr_final), "genes ×", ncol(expr_final), "samples\n")
cat("保存到 Data/GSE208668_processed.RData\n")
