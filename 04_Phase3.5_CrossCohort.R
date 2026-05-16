# ============================================================
# 04_Phase3.5_CrossCohort.R
# ============================================================
# 对应 Methods §2.5 Cross-cohort directional consistency check
#
# 输入：
#   Data/GSE208668_processed.RData (基准队列)
#   在线下载 GSE39445 / GSE98566 / GSE56931
# 输出：
#   Phase3.5_CrossCohort/GSE39445_target_genes.csv
#   Phase3.5_CrossCohort/GSE98566_target_genes.csv
#   Phase3.5_CrossCohort/GSE56931_target_genes.csv
#   Phase3.5_CrossCohort/cross_cohort_summary.csv  (Table S7)
#
# 注意：3 个验证队列均为急性睡眠剥夺（不同病理生理），
#       只能做 directional consistency check，不能叫 "validation"。
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")
dir.create("Phase3.5_CrossCohort/raw", recursive = TRUE, showWarnings = FALSE)

library(GEOquery)
library(limma)
library(data.table)

# ============================================================
# Step 0：在 GSE208668 中先获取 9 个目标基因的基准 logFC
# ============================================================
load("Data/GSE208668_processed.RData")

target_genes <- c("SPTLC1", "PTGES3", "HLA-G",       # ML Tier 1+2
                  "NFE2L2", "SOD2", "HMOX1",          # OxStress
                  "TNFAIP3", "P2RX7", "FKBP5")        # Hypothesis-driven

design_208668 <- model.matrix(~ 0 + group)
colnames(design_208668) <- c("control", "insomnia")
fit_208668 <- lmFit(expr_final, design_208668)
contrast_208668 <- makeContrasts(insomnia - control, levels = design_208668)
fit2_208668 <- eBayes(contrasts.fit(fit_208668, contrast_208668))
tt_208668 <- topTable(fit2_208668, number = Inf)

gse208668_results <- data.frame()
for (gene in target_genes) {
  if (gene %in% rownames(tt_208668)) {
    row <- tt_208668[gene, ]
    gse208668_results <- rbind(gse208668_results, data.frame(
      Gene = gene,
      GSE208668_logFC     = round(row$logFC, 4),
      GSE208668_P         = signif(row$P.Value, 3),
      GSE208668_Direction = ifelse(row$logFC > 0, "UP", "DOWN"),
      stringsAsFactors = FALSE
    ))
  }
}
cat("=== GSE208668 基准（9 基因均为 UP）===\n")
print(gse208668_results)

# ============================================================
# Cohort 1：GSE39445（急性睡眠限制, n=26, Möller-Levet 2013）
# ============================================================
cat("\n========== GSE39445 ==========\n")
options(timeout = 3600)

gse39445 <- getGEO("GSE39445", GSEMatrix = TRUE, getGPL = TRUE,
                   destdir = "Phase3.5_CrossCohort/raw")
eset <- gse39445[[1]]
pheno39 <- pData(eset)

# 提取睡眠条件 + 被试 ID
sleep_protocol <- pheno39$`sleepprotocol:ch1`
subject_id     <- gsub("subject: ", "", pheno39$`subject:ch1`)

# 排除 NA 被试
keep <- !is.na(subject_id) & subject_id != "NA" & subject_id != ""
expr39      <- exprs(eset)[, keep]
sleep_cond  <- factor(sleep_protocol[keep],
                      levels = c("Sleep Extension", "Sleep Restriction"))
subj        <- factor(subject_id[keep])

cat("有效样本:", sum(keep), "\n")
cat("条件分布:\n");  print(table(sleep_cond))

# limma + duplicateCorrelation（同一受试者多次采样需校正）
design39 <- model.matrix(~ 0 + sleep_cond)
colnames(design39) <- c("Extension", "Restriction")

cat("计算 duplicateCorrelation ...\n")
corfit <- duplicateCorrelation(expr39, design39, block = subj)
cat("组内相关系数:", corfit$consensus.correlation, "\n")

fit39  <- lmFit(expr39, design39, block = subj,
                correlation = corfit$consensus.correlation)
contrast39 <- makeContrasts(Restriction - Extension, levels = design39)
fit39_2 <- eBayes(contrasts.fit(fit39, contrast39))

# 提取目标基因
feat3      <- fData(eset)
gene_symbol <- feat3$GENE_SYMBOL
tt39_all   <- topTable(fit39_2, number = Inf)
tt39_all$GENE_SYMBOL <- gene_symbol[match(rownames(tt39_all), rownames(expr39))]

results_39445 <- data.frame()
for (gene in target_genes) {
  rows <- tt39_all[which(tt39_all$GENE_SYMBOL == gene), ]
  if (nrow(rows) > 0) {
    best <- rows[which.min(rows$P.Value), ]
    results_39445 <- rbind(results_39445, data.frame(
      Gene = gene,
      GSE39445_logFC     = round(best$logFC, 4),
      GSE39445_P         = signif(best$P.Value, 3),
      GSE39445_adj_P     = signif(best$adj.P.Val, 3),
      GSE39445_Direction = ifelse(best$logFC > 0, "UP", "DOWN"),
      N_probes = nrow(rows),
      stringsAsFactors = FALSE
    ))
  } else {
    results_39445 <- rbind(results_39445, data.frame(
      Gene = gene, GSE39445_logFC = NA, GSE39445_P = NA,
      GSE39445_adj_P = NA, GSE39445_Direction = "NOT_FOUND", N_probes = 0
    ))
  }
}
print(results_39445)
write.csv(results_39445, "Phase3.5_CrossCohort/GSE39445_target_genes.csv",
          row.names = FALSE)

# ============================================================
# Cohort 2：GSE98566（62h 总睡眠剥夺, n=17, Uyhelji 2018）
# ============================================================
cat("\n========== GSE98566 ==========\n")
gse98566 <- getGEO("GSE98566", GSEMatrix = TRUE, getGPL = TRUE,
                   destdir = "Phase3.5_CrossCohort/raw")
eset2 <- gse98566[[1]]
pheno98 <- pData(eset2)

subj_group98 <- pheno98$`subject group:ch1`
subj_id98    <- gsub("subject id: ", "", pheno98$`subject id:ch1`)

# 基因符号在 gene_assignment 列中（格式：NM_xxx // GENE_SYMBOL // ...）
feat2 <- fData(eset2)
extract_gene <- function(x) {
  if (is.na(x) || x == "" || x == "---") return(NA)
  parts <- strsplit(x, " // ")[[1]]
  if (length(parts) >= 2) return(trimws(parts[2]))
  return(NA)
}
gene_symbol98 <- sapply(feat2$gene_assignment, extract_gene, USE.NAMES = FALSE)

keep98 <- !is.na(subj_group98) & !is.na(subj_id98)
expr98 <- exprs(eset2)[, keep98]
grp98  <- factor(subj_group98[keep98], levels = c("Control", "Sleep Deprived"))
subj98 <- factor(subj_id98[keep98])

design98 <- model.matrix(~ 0 + grp98)
colnames(design98) <- c("Control", "SleepDeprived")

corfit98 <- duplicateCorrelation(expr98, design98, block = subj98)
fit98    <- lmFit(expr98, design98, block = subj98,
                  correlation = corfit98$consensus.correlation)
contrast98 <- makeContrasts(SleepDeprived - Control, levels = design98)
fit98_2  <- eBayes(contrasts.fit(fit98, contrast98))

tt98_all <- topTable(fit98_2, number = Inf)
tt98_all$GENE_SYMBOL <- gene_symbol98[match(rownames(tt98_all), rownames(expr98))]

results_98566 <- data.frame()
for (gene in target_genes) {
  rows <- tt98_all[which(tt98_all$GENE_SYMBOL == gene), ]
  if (nrow(rows) > 0) {
    best <- rows[which.min(rows$P.Value), ]
    results_98566 <- rbind(results_98566, data.frame(
      Gene = gene,
      GSE98566_logFC     = round(best$logFC, 4),
      GSE98566_P         = signif(best$P.Value, 3),
      GSE98566_adj_P     = signif(best$adj.P.Val, 3),
      GSE98566_Direction = ifelse(best$logFC > 0, "UP", "DOWN"),
      N_probes = nrow(rows),
      stringsAsFactors = FALSE
    ))
  } else {
    results_98566 <- rbind(results_98566, data.frame(
      Gene = gene, GSE98566_logFC = NA, GSE98566_P = NA,
      GSE98566_adj_P = NA, GSE98566_Direction = "NOT_FOUND", N_probes = 0
    ))
  }
}
print(results_98566)
write.csv(results_98566, "Phase3.5_CrossCohort/GSE98566_target_genes.csv",
          row.names = FALSE)

# ============================================================
# Cohort 3：GSE56931（38h 持续清醒, n=14, Pellegrino 2014）
# ============================================================
cat("\n========== GSE56931 ==========\n")
# 平台 GPL10379 注释需手动下载（GEOquery getGEO 易失败）
# 这里直接读取已下载的 series_matrix
gse56931 <- getGEO(filename = "Phase3.5_CrossCohort/raw/GSE56931_series_matrix.txt.gz",
                   getGPL = FALSE)
pheno56 <- pData(gse56931)

bio_group  <- gsub("biological group: ", "", pheno56$`biological group:ch1`)
patient_id <- gsub("patient: ", "", pheno56$`patient:ch1`)

# 加载 GPL10379 平台注释（手动从 NCBI 获取）
gpl_annot <- fread("Phase3.5_CrossCohort/raw/GPL10379_annot.txt", header = TRUE)

# 只保留 baseline + sleep deprivation（排除 recovery 期）
keep56 <- bio_group %in% c("baseline", "sleep deprivation")
expr56 <- exprs(gse56931)[, keep56]
grp56  <- factor(bio_group[keep56], levels = c("baseline", "sleep deprivation"))
subj56 <- factor(patient_id[keep56])

design56 <- model.matrix(~ 0 + grp56)
colnames(design56) <- c("Baseline", "SleepDep")

corfit56 <- duplicateCorrelation(expr56, design56, block = subj56)
fit56    <- lmFit(expr56, design56, block = subj56,
                  correlation = corfit56$consensus.correlation)
contrast56 <- makeContrasts(SleepDep - Baseline, levels = design56)
fit56_2  <- eBayes(contrasts.fit(fit56, contrast56))

tt56_all <- topTable(fit56_2, number = Inf)
gene_map <- gpl_annot[, c("ID", "GeneSymbol")]
gene_map$ID <- as.character(gene_map$ID)
tt56_all$GENE_SYMBOL <- gene_map$GeneSymbol[match(rownames(tt56_all), gene_map$ID)]

results_56931 <- data.frame()
for (gene in target_genes) {
  rows <- tt56_all[which(tt56_all$GENE_SYMBOL == gene), ]
  if (nrow(rows) > 0) {
    best <- rows[which.min(rows$P.Value), ]
    results_56931 <- rbind(results_56931, data.frame(
      Gene = gene,
      GSE56931_logFC     = round(best$logFC, 4),
      GSE56931_P         = signif(best$P.Value, 3),
      GSE56931_adj_P     = signif(best$adj.P.Val, 3),
      GSE56931_Direction = ifelse(best$logFC > 0, "UP", "DOWN"),
      N_probes = nrow(rows),
      stringsAsFactors = FALSE
    ))
  } else {
    results_56931 <- rbind(results_56931, data.frame(
      Gene = gene, GSE56931_logFC = NA, GSE56931_P = NA,
      GSE56931_adj_P = NA, GSE56931_Direction = "NOT_FOUND", N_probes = 0
    ))
  }
}
print(results_56931)
write.csv(results_56931, "Phase3.5_CrossCohort/GSE56931_target_genes.csv",
          row.names = FALSE)

# ============================================================
# 汇总：4 数据集 × 9 基因
# ============================================================
cat("\n========== 跨队列汇总 ==========\n")
summary_table <- gse208668_results[, c("Gene", "GSE208668_logFC", "GSE208668_Direction")]
summary_table <- merge(summary_table,
                       results_39445[, c("Gene", "GSE39445_logFC", "GSE39445_P", "GSE39445_Direction")],
                       by = "Gene", all.x = TRUE)
summary_table <- merge(summary_table,
                       results_98566[, c("Gene", "GSE98566_logFC", "GSE98566_P", "GSE98566_Direction")],
                       by = "Gene", all.x = TRUE)
summary_table <- merge(summary_table,
                       results_56931[, c("Gene", "GSE56931_logFC", "GSE56931_P", "GSE56931_Direction")],
                       by = "Gene", all.x = TRUE)

# 方向一致性统计
summary_table$Consistent_N <- apply(summary_table, 1, function(row) {
  ref_dir <- row["GSE208668_Direction"]
  n_consistent <- 0;  n_available <- 0
  for (col in c("GSE39445_Direction", "GSE98566_Direction", "GSE56931_Direction")) {
    if (!is.na(row[col]) & row[col] != "NOT_FOUND") {
      n_available <- n_available + 1
      if (row[col] == ref_dir) n_consistent <- n_consistent + 1
    }
  }
  paste0(n_consistent, "/", n_available)
})

# 添加 Tier 标签
summary_table$Tier <- ifelse(summary_table$Gene == "SPTLC1", "1_ML_3algo",
                      ifelse(summary_table$Gene %in% c("PTGES3", "HLA-G"), "2_ML_2algo",
                      ifelse(summary_table$Gene %in% c("NFE2L2", "SOD2", "HMOX1"), "3_OxStress",
                             "4_Hypothesis")))

summary_table <- summary_table[order(summary_table$Tier, summary_table$Gene), ]

cat("\n=== Table S7：跨队列方向一致性 ===\n")
print(summary_table, row.names = FALSE)

write.csv(summary_table, "Phase3.5_CrossCohort/cross_cohort_summary.csv",
          row.names = FALSE)

cat("\n========== Phase 3.5 完成 ==========\n")
cat("结论：慢性失眠中 9 基因均上调 (logFC +0.94 to +2.81)，\n")
cat("      但在 3 个急性睡眠剥夺队列中 logFC 接近 0（-0.22 to +0.10），\n")
cat("      方向多不一致 → 提示慢性 vs 急性是不同的病理生理过程。\n")
cat("⚠️  论文中只能称 'directional consistency check'，不能叫 'validation'。\n")
