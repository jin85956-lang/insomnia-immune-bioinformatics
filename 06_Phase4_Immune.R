# ============================================================
# 06_Phase4_Immune.R
# ============================================================
# 对应 Methods §2.7 Immune microenvironment analysis (CIBERSORTx)
#
# 输入：
#   Data/GSE208668_processed.RData
#   Phase1_GSVA/GSVA_results.RData
# 输出（两阶段）：
#   阶段 A — CIBERSORTx 输入准备
#     Phase4_Immune/CIBERSORTx_input.txt  (上传到 CIBERSORTx 网站)
#   阶段 B — CIBERSORTx 结果分析（CIBERSORTx 网站跑完后再继续）
#     Phase4_Immune/CIBERSORTx_Job460_Results.csv  (网站返回的 CSV，手动放到本地)
#     Phase4_Immune/immune_cell_diff.csv           (22 细胞差异统计, Table S8)
#     Phase4_Immune/GSVA_immune_correlation.csv    (5 通路 × 10 细胞 Spearman)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")

library(ggplot2)
library(reshape2)
library(pheatmap)

# ============================================================
# 阶段 A — 准备 CIBERSORTx 输入
# ============================================================
load("Data/GSE208668_processed.RData")

# CIBERSORTx 要求：基因符号 + 线性 expression（不是 log2）
# GSE208668 的 expr_final 是 log2 微阵列数据，需要反 log
expr_linear <- 2^expr_final
cat("反 log2 后表达值范围:", range(expr_linear), "\n")
cat("中位数:", round(median(expr_linear), 1), "\n")

# 第一列必须是 GeneSymbol，其余列是样本（tab 分隔）
cibersort_input <- data.frame(
  GeneSymbol = rownames(expr_linear),
  as.data.frame(expr_linear),
  check.names = FALSE
)

write.table(cibersort_input,
            file = "Phase4_Immune/CIBERSORTx_input.txt",
            sep  = "\t",
            quote = FALSE, row.names = FALSE)

cat("\n========== CIBERSORTx 输入文件已生成 ==========\n")
cat("文件路径: Phase4_Immune/CIBERSORTx_input.txt\n")
cat("基因数:", nrow(cibersort_input), "\n")
cat("样本数:", ncol(cibersort_input) - 1, "\n")
cat("文件大小:", round(file.size("Phase4_Immune/CIBERSORTx_input.txt")/1e6, 1), "MB\n")
cat("\n下一步：\n")
cat("  1. 登录 https://cibersortx.stanford.edu/\n")
cat("  2. 选择 Impute Cell Fractions 模块\n")
cat("  3. Signature matrix: LM22\n")
cat("  4. Mixture file: 上传 CIBERSORTx_input.txt\n")
cat("  5. Permutations: 100\n")
cat("  6. 跑完后下载 Job Results CSV，重命名为 CIBERSORTx_Job460_Results.csv\n")
cat("     存放到 Phase4_Immune/ 目录下，继续运行下面的阶段 B 代码\n\n")

# ============================================================
# 阶段 B — CIBERSORTx 结果分析
# ============================================================
# 以下代码在 CIBERSORTx 跑完、CSV 下载回本地后再执行
# 如果 CSV 还未生成，下面会报错——这是正常的
# ============================================================

cibersort_file <- "Phase4_Immune/CIBERSORTx_Job460_Results.csv"
if (!file.exists(cibersort_file)) {
  cat("⚠️  CIBERSORTx 结果文件还不存在。\n")
  cat("    请先到 CIBERSORTx 网站跑完，下载 CSV 后再继续。\n")
  stop("等待 CIBERSORTx 结果")
}

# ---- 1. 加载 CIBERSORTx 输出 ----
ciber <- read.csv(cibersort_file, row.names = 1, check.names = FALSE)
cat("CIBERSORTx 结果维度:", dim(ciber), "\n")
cat("列名:\n");  print(colnames(ciber))
# 22 细胞 + P-value + Correlation + RMSE = 25 列

# ---- 2. 分离细胞比例与 QC 信息 ----
cell_fractions <- ciber[, 1:22]
qc_info        <- ciber[, 23:25]

cat("\n=== CIBERSORTx QC ===\n")
cat("P-value < 0.05 的样本数:", sum(qc_info$`P-value` < 0.05),
    "/", nrow(qc_info), "\n")
cat("P-value 范围:", range(qc_info$`P-value`), "\n")
cat("Correlation 范围:", range(qc_info$Correlation), "\n")
cat("RMSE 范围:", range(qc_info$RMSE), "\n")

# ---- 3. 添加分组信息 ----
sample_ids   <- rownames(cell_fractions)
group_match  <- group[match(sample_ids, colnames(expr_final))]
cat("\n分组匹配:");  print(table(group_match))

# ---- 4. 22 种细胞 × 失眠 vs 对照差异分析 ----
cat("\n========================================\n")
cat("  22 种免疫细胞差异分析\n")
cat("========================================\n\n")

results_immune <- data.frame(
  CellType      = character(),
  Mean_Insomnia = numeric(),
  Mean_Control  = numeric(),
  Diff          = numeric(),
  Direction     = character(),
  Wilcox_P      = numeric(),
  Cohen_d       = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:22) {
  cell_name <- colnames(cell_fractions)[i]
  vals_ins  <- cell_fractions[group_match == "yes", i]
  vals_ctrl <- cell_fractions[group_match == "no",  i]
  
  wt <- wilcox.test(vals_ins, vals_ctrl)
  
  pooled_sd <- sqrt(((length(vals_ins)-1)*sd(vals_ins)^2 + 
                     (length(vals_ctrl)-1)*sd(vals_ctrl)^2) /
                    (length(vals_ins) + length(vals_ctrl) - 2))
  cohens_d  <- ifelse(pooled_sd == 0, 0,
                      (mean(vals_ins) - mean(vals_ctrl)) / pooled_sd)
  
  results_immune <- rbind(results_immune, data.frame(
    CellType      = cell_name,
    Mean_Insomnia = round(mean(vals_ins),  4),
    Mean_Control  = round(mean(vals_ctrl), 4),
    Diff          = round(mean(vals_ins) - mean(vals_ctrl), 4),
    Direction     = ifelse(mean(vals_ins) > mean(vals_ctrl), "失眠↑", "失眠↓"),
    Wilcox_P      = round(wt$p.value, 4),
    Cohen_d       = round(cohens_d, 3),
    stringsAsFactors = FALSE
  ))
}

results_immune <- results_immune[order(results_immune$Wilcox_P), ]
results_immune$Sig <- ifelse(results_immune$Wilcox_P < 0.01, "**",
                      ifelse(results_immune$Wilcox_P < 0.05, "*",
                      ifelse(results_immune$Wilcox_P < 0.1,  ".",  "")))

cat("=== 完整结果（按 p 值排序）===\n")
print(results_immune, row.names = FALSE)

sig_cells <- results_immune[results_immune$Wilcox_P < 0.05, ]
cat("\n=== 显著差异细胞（p<0.05），共", nrow(sig_cells), "种 ===\n")
print(sig_cells, row.names = FALSE)

write.csv(results_immune, "Phase4_Immune/immune_cell_diff.csv", row.names = FALSE)

# ---- 5. GSVA 通路 × 免疫细胞 Spearman 相关 ----
cat("\n========================================\n")
cat("  GSVA 通路 vs 免疫细胞相关性\n")
cat("========================================\n\n")

load("Phase1_GSVA/GSVA_results.RData")

sig_cells_names <- results_immune$CellType[results_immune$Wilcox_P < 0.05]
sig_immune <- cell_fractions[, sig_cells_names]

# GSVA 分数对齐到 CIBERSORTx 样本
gsva_mat <- gsva_scores[, match(rownames(cell_fractions), colnames(gsva_scores))]

cor_results <- data.frame()
cor_matrix <- matrix(NA, nrow = nrow(gsva_mat), ncol = ncol(sig_immune))
p_matrix   <- cor_matrix
rownames(cor_matrix) <- rownames(gsva_mat)
colnames(cor_matrix) <- colnames(sig_immune)

for (i in 1:nrow(gsva_mat)) {
  for (j in 1:ncol(sig_immune)) {
    ct <- cor.test(as.numeric(gsva_mat[i, ]), sig_immune[, j],
                   method = "spearman")
    cor_matrix[i, j] <- ct$estimate
    p_matrix[i, j]   <- ct$p.value
    cor_results <- rbind(cor_results, data.frame(
      Pathway  = rownames(gsva_mat)[i],
      CellType = colnames(sig_immune)[j],
      Rho      = round(ct$estimate, 3),
      P_value  = round(ct$p.value,  4),
      stringsAsFactors = FALSE
    ))
  }
}

sig_cor <- cor_results[cor_results$P_value < 0.05, ]
sig_cor <- sig_cor[order(-abs(sig_cor$Rho)), ]
cat("显著相关对数（p<0.05）:", nrow(sig_cor), "\n")
print(head(sig_cor, 15))

write.csv(cor_results, "Phase4_Immune/GSVA_immune_correlation.csv", row.names = FALSE)

cat("\n========== Phase 4 免疫分析完成 ==========\n")
cat("⚠️  关键警告：\n")
cat("  1. CIBERSORTx QC 失败的 8 个样本全部来自 Control 组（潜在偏倚）\n")
cat("  2. M1 macrophages 和 DC resting 在失眠组均值=0（floor effect）\n")
cat("  3. 论文中只能称 'relative estimates'，不能说 'absolute cell counts'\n")
