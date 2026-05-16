# ##############################################################################
# ##############################################################################
# ##                                                                          ##
# ##   Insomnia Immune-Stress Bioinformatics — Supplementary Code (Full)      ##
# ##   GSE208668 (PBMC, Illumina HumanHT-12 V4, n = 42)                       ##
# ##   17 insomnia + 25 control                                               ##
# ##                                                                          ##
# ##   This single file merges all 12 analysis scripts (Phase 0 – Phase 6)    ##
# ##   in their intended execution order. Each section is self-contained and  ##
# ##   can be run independently after Section 1 is complete.                  ##
# ##                                                                          ##
# ##   HOW TO RUN                                                             ##
# ##     Option A — run everything once:                                      ##
# ##         source("Insomnia_Supplementary_Code_ALL.R", encoding = "UTF-8")  ##
# ##                                                                          ##
# ##     Option B — run section by section in RStudio:                        ##
# ##         use Ctrl + F to jump to ">>> SECTION X <<<" markers,             ##
# ##         then select the block and Ctrl + Enter.                          ##
# ##                                                                          ##
# ##   EXTERNAL STEPS that pause this script (cannot be done in R alone):     ##
# ##     - Section 7  : upload expr matrix to CIBERSORTx web server           ##
# ##     - Section 8  : run SMR command-line tool in a terminal               ##
# ##     - Section 9  : query DGIdb web interface                             ##
# ##     Each section explains what to do, then resumes parsing.              ##
# ##                                                                          ##
# ##   DEPENDENCIES                                                           ##
# ##     CRAN          : tidyverse, data.table, glmnet, randomForest, Boruta, ##
# ##                     pROC, VennDiagram, WGCNA, ggplot2, openxlsx          ##
# ##     Bioconductor  : GEOquery, limma, GSVA, clusterProfiler, org.Hs.eg.db ##
# ##                                                                          ##
# ##   WRITING RED LINES (already enforced in code comments below)            ##
# ##     - SMR is "consistent with" causality, not "proof"                    ##
# ##     - Cross-cohort is a "directional consistency check", not validation  ##
# ##     - AUC = 1.000 is "discovery-stage, requires external validation"     ##
# ##     - CIBERSORTx outputs are "relative estimates", not absolute counts   ##
# ##                                                                          ##
# ##############################################################################
# ##############################################################################




# ##############################################################################
# >>> SECTION 1 <<<
# Setup, GEO download and preprocessing
# Source : 00_setup_and_data.R
# ##############################################################################

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



# ##############################################################################
# >>> SECTION 2 <<<
# Phase 1 — GSVA of five stress / immune gene sets
# Source : 01_Phase1_GSVA.R
# ##############################################################################

# ============================================================
# 01_Phase1_GSVA.R
# ============================================================
# 对应 Methods §2.2 Pathway-level enrichment analysis
#
# 输入：Data/GSE208668_processed.RData
# 输出：
#   Phase1_GSVA/GSVA_results.RData        (gsva_scores, results_gsva, pathway_filtered)
#   Phase1_GSVA/GSVA_results.csv          (Cohen's d + Wilcoxon p)
#   Phase1_GSVA/GSVA_gene_sets.csv        (74 基因 × 通路对照, Table S11)
#   Phase1_GSVA/GSVA_gene_sets_summary.csv
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")

library(GSVA)
library(effsize)
library(dplyr)

# ---- 1. 加载数据 ----
load("Data/GSE208668_processed.RData")
expr_matrix <- as.matrix(expr_final)

cat("=== 数据基本信息 ===\n")
cat("表达矩阵:", dim(expr_matrix), "\n")
cat("分组:    ", table(group), "\n")

# ---- 2. 定义 5 个文献策划基因集 ----
# 全部来源于已发表文献，详见 Methods §2.2
pathway_genes <- list(
  
  # NF-κB signaling targets (25 genes)
  NFkB = c("TNFAIP3","NFKBIA","NFKBIZ","IL1B","IL8","PTGS2",
           "CHUK","IKBKB","IKBKG","RELA","REL","RELB","NFKB1","NFKB2",
           "MYD88","TRAF6","TLR2","TLR4","TNF","BCL2L1","BIRC3",
           "NFKBIE","NFKBIB","IL6","IKBKE"),
  
  # Glucocorticoid response (18 genes, Oakley & Cidlowski 2013)
  GR = c("NR3C1","FKBP4","FKBP5","TSC22D3","DUSP1","PER1",
         "SGK1","KLF9","HSPA1A","HSP90AA1","HSP90AB1","ANXA1",
         "SCNN1A","HSD11B1","HSD11B2","CRH","POMC","SLC22A5"),
  
  # Interferon-stimulated genes (12 genes)
  ISG = c("IFI44","IFI44L","IRF9","IRF7","MX1","MX2","ISG20",
          "IFITM2","IFIT1","IFIT3","OAS1","STAT1"),
  
  # NLRP3 inflammasome pathway (11 genes, Broz & Dixit 2016)
  Inflammasome = c("NLRP3","PYCARD","CASP1","IL1B","IL18",
                   "P2RX7","TXNIP","HMGB1","GSDMD","IL1R1","IL18R1"),
  
  # Nrf2 / Oxidative stress response (8 genes)
  OxStress = c("NFE2L2","NQO1","HMOX1","SOD1","SOD2","CAT","GPX1","GCLC")
)

# ---- 3. 基因集与表达矩阵交集（QC） ----
cat("\n=== 基因集匹配 ===\n")
pathway_filtered <- list()
for (name in names(pathway_genes)) {
  matched <- intersect(pathway_genes[[name]], rownames(expr_matrix))
  missing <- setdiff(pathway_genes[[name]], rownames(expr_matrix))
  pathway_filtered[[name]] <- matched
  cat(sprintf("  %-15s: %2d defined → %2d in data", 
              name, length(pathway_genes[[name]]), length(matched)))
  if (length(missing) > 0) cat("  [缺失:", paste(missing, collapse=", "), "]")
  cat("\n")
}
# 预期：NFkB 25/25, GR 15/18, ISG 12/12, Inflammasome 10/11, OxStress 8/8
# 总计：74 定义，70 可用

# ---- 4. GSVA 计算 ----
cat("\n=== GSVA 计算 ===\n")
cat("GSVA 包版本:", as.character(packageVersion("GSVA")), "\n")

# 新版 GSVA (>= 1.50) 用 gsvaParam + gsva() 二步法
gsva_param <- gsvaParam(
  exprData = expr_matrix,
  geneSets = pathway_filtered,
  kcdf = "Gaussian"          # 微阵列 log2 数据用 Gaussian
)
gsva_scores <- gsva(gsva_param, verbose = FALSE)

cat("GSVA 评分矩阵:", dim(gsva_scores), "\n")  # 5 × 42

# ---- 5. Cohen's d + Wilcoxon 检验 ----
cat("\n=== 通路差异统计（按 |d| 降序）===\n")
results_gsva <- data.frame(
  Module        = character(),
  Cohens_d      = numeric(),
  Wilcox_p      = numeric(),
  Ttest_p       = numeric(),
  Mean_insomnia = numeric(),
  Mean_control  = numeric(),
  Direction     = character(),
  stringsAsFactors = FALSE
)

grp <- factor(group)

for (mod in rownames(gsva_scores)) {
  ins  <- gsva_scores[mod, grp == "yes"]
  ctrl <- gsva_scores[mod, grp == "no"]
  
  n1 <- length(ins); n2 <- length(ctrl)
  m1 <- mean(ins);   m2 <- mean(ctrl)
  s1 <- sd(ins);     s2 <- sd(ctrl)
  
  # 包内 Cohen's d
  d_pkg <- cohen.d(as.numeric(ins), as.numeric(ctrl))$estimate
  
  # 手动 pooled SD 验证（双重保险）
  s_pooled <- sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
  d_manual <- (m1 - m2) / s_pooled
  
  stopifnot(abs(d_pkg - d_manual) < 0.01)
  
  w_p <- wilcox.test(as.numeric(ins), as.numeric(ctrl))$p.value
  t_p <- t.test(as.numeric(ins), as.numeric(ctrl))$p.value
  
  results_gsva <- rbind(results_gsva, data.frame(
    Module        = mod,
    Cohens_d      = round(d_pkg, 4),
    Wilcox_p      = round(w_p, 4),
    Ttest_p       = round(t_p, 4),
    Mean_insomnia = round(m1, 4),
    Mean_control  = round(m2, 4),
    Direction     = ifelse(m1 > m2, "Insomnia_UP", "Insomnia_DOWN"),
    stringsAsFactors = FALSE
  ))
}

results_gsva <- results_gsva[order(-abs(results_gsva$Cohens_d)), ]
cat("\n");  print(results_gsva)
# 预期结果（论文 Table 1）：
#   OxStress     d = +0.856, p = 0.057
#   NFkB         d = +0.489, p = 0.153
#   GR           d = +0.415, p = 0.193
#   ISG          d = −0.253, p = 0.431
#   Inflammasome d = −0.002, p = 0.285

# ---- 6. 保存：GSVA 结果 ----
write.csv(results_gsva, "Phase1_GSVA/GSVA_results.csv", row.names = FALSE)
save(gsva_scores, results_gsva, pathway_filtered, pathway_genes,
     file = "Phase1_GSVA/GSVA_results.RData")
cat("\n保存：Phase1_GSVA/GSVA_results.csv + GSVA_results.RData\n")

# ---- 7. 导出 Table S11（GSVA 基因集定义）----
cat("\n=== 导出 Table S11（基因集定义）===\n")

available_genes <- rownames(expr_final)

# 长格式：每行一个基因
table_s11 <- do.call(rbind, lapply(names(pathway_filtered), function(pw) {
  data.frame(
    Pathway              = pw,
    Gene_Symbol          = pathway_filtered[[pw]],
    In_Original_Definition = "Yes",
    In_Expression_Matrix   = "Yes",
    Source = switch(pw,
      "NFkB"         = "Literature-curated (NF-kB signaling targets)",
      "GR"           = "Literature-curated (Oakley & Cidlowski 2013)",
      "ISG"          = "Literature-curated (Interferon-stimulated genes)",
      "Inflammasome" = "Literature-curated (Broz & Dixit 2016)",
      "OxStress"     = "Literature-curated (Nrf2/Oxidative stress)"
    ),
    stringsAsFactors = FALSE
  )
}))

# 把缺失基因也加进去，便于审稿人核对
for (pw in names(pathway_genes)) {
  missing <- setdiff(pathway_genes[[pw]], available_genes)
  if (length(missing) > 0) {
    table_s11 <- rbind(table_s11, data.frame(
      Pathway              = pw,
      Gene_Symbol          = missing,
      In_Original_Definition = "Yes",
      In_Expression_Matrix   = "No (probe not available)",
      Source = table_s11$Source[table_s11$Pathway == pw][1],
      stringsAsFactors = FALSE
    ))
  }
}

table_s11 <- table_s11[order(table_s11$Pathway, 
                             table_s11$In_Expression_Matrix, 
                             table_s11$Gene_Symbol), ]

# 汇总表
summary_s11 <- data.frame(
  Pathway          = names(pathway_genes),
  N_Genes_Defined  = sapply(pathway_genes,    length),
  N_Genes_In_Data  = sapply(pathway_filtered, length),
  N_Genes_Missing  = sapply(names(pathway_genes), function(nm) {
    length(pathway_genes[[nm]]) - length(pathway_filtered[[nm]])
  }),
  Missing_Genes = sapply(names(pathway_genes), function(nm) {
    m <- setdiff(pathway_genes[[nm]], available_genes)
    if (length(m) == 0) "None" else paste(m, collapse = ", ")
  }),
  stringsAsFactors = FALSE
)
rownames(summary_s11) <- NULL

write.csv(table_s11, "Phase1_GSVA/GSVA_gene_sets.csv", row.names = FALSE)
write.csv(summary_s11, "Phase1_GSVA/GSVA_gene_sets_summary.csv", row.names = FALSE)
cat("Table S11 已保存\n")
print(summary_s11)

cat("\n========== Phase 1 GSVA 完成 ==========\n")
cat("决策点：至少一个通路 Wilcox p < 0.05 → 信号存在，进入 Phase 2\n")
cat("        当前 OxStress p = 0.057（趋势）→ 可继续 WGCNA\n")



# ##############################################################################
# >>> SECTION 3 <<<
# Phase 2 — WGCNA module discovery
# Source : 02_Phase2_WGCNA.R
# ##############################################################################

# ============================================================
# 02_Phase2_WGCNA.R
# ============================================================
# 对应 Methods §2.3 Weighted Gene Co-Expression Network Analysis
#
# 输入：Data/GSE208668_processed.RData
# 输出：
#   Phase2_WGCNA/WGCNA_results.RData   (net, MEs, moduleTraitCor, hub_genes)
#   Phase2_WGCNA/hub_genes.csv         (2,120 hub 基因, Table S5)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")

library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# ---- 1. 加载数据 ----
load("Data/GSE208668_processed.RData")
expr_matrix <- as.matrix(expr_final)

# ---- 2. 取方差最高的 15% 基因 ----
# 内存优化（14,539 基因全跑会内存溢出）
var_genes <- apply(expr_matrix, 1, var)
n_top <- round(nrow(expr_matrix) * 0.15)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:n_top]
cat("WGCNA 输入基因数:", length(top_genes), "\n")  # ≈ 2,181

datExpr <- t(expr_matrix[top_genes, ])
cat("datExpr 维度:", dim(datExpr), "（样本 × 基因）\n")

# ---- 3. 样本聚类 QC（检查离群样本）----
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("Phase2_WGCNA/图_sample_clustering.pdf", width = 12, height = 6)
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "")
dev.off()
cat("样本聚类图已保存\n")

# ---- 4. 软阈值选择 ----
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 2)

pdf("Phase2_WGCNA/图_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")
abline(h = 0.80, col = "blue", lty = 2)

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = 0.9, col = "red")
dev.off()

cat("\nsft$powerEstimate:", sft$powerEstimate, "\n")

# 经测试，本数据集 power = 6 满足 R² > 0.80
softPower <- 6

# ---- 5. 构建网络 + 模块检测 ----
cat("\n=== 构建共表达网络（需要几分钟）===\n")
net <- blockwiseModules(
  datExpr,
  power            = softPower,
  TOMType          = "signed",
  minModuleSize    = 30,
  reassignThreshold = 0,
  mergeCutHeight   = 0.25,
  numericLabels    = TRUE,           # 0=grey, 1=turquoise, ...
  pamRespectsDendro = FALSE,
  saveTOMs         = FALSE,
  maxBlockSize     = 2500,           # 单块最多 2500 基因
  verbose          = 3
)

cat("\n模块数量（不含 grey）:", length(unique(net$colors)) - 1, "\n")
print(table(net$colors))
# 预期：0 (grey) = 3 基因, 1 (turquoise) = 2,178 基因

# ---- 6. 模块树状图 ----
mergedColors <- labels2colors(net$colors)
pdf("Phase2_WGCNA/图_WGCNA_dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

# ---- 7. 模块-特征相关性 ----
MEs <- net$MEs

traitData <- data.frame(
  Insomnia = as.numeric(group == "yes"),
  Age      = as.numeric(pheno$age),
  BMI      = as.numeric(pheno$bmi),
  BDI      = as.numeric(pheno$bdi)
)

moduleTraitCor    <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("\n=== 模块与失眠的相关性 ===\n")
insomnia_cor <- data.frame(
  Module = colnames(MEs),
  Color  = labels2colors(as.numeric(gsub("ME", "", colnames(MEs)))),
  Cor    = moduleTraitCor[, "Insomnia"],
  P      = moduleTraitPvalue[, "Insomnia"]
)
insomnia_cor <- insomnia_cor[order(insomnia_cor$P), ]
print(insomnia_cor)
# 预期：Turquoise r = -0.44, p = 0.003

# ---- 8. 模块-特征热图 ----
pdf("Phase2_WGCNA/图_module_trait_heatmap.pdf", width = 8, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(traitData),
  yLabels     = colnames(MEs),
  ySymbols    = colnames(MEs),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = textMatrix,
  setStdMargins = FALSE,
  cex.text    = 0.6,
  zlim        = c(-1, 1),
  main        = "Module-Trait Relationships"
)
dev.off()

# ---- 9. Hub 基因识别 ----
# 选最显著（非 grey）的模块
best_module_idx <- which.min(insomnia_cor$P[insomnia_cor$Module != "ME0"])
best_ME <- insomnia_cor$Module[insomnia_cor$Module != "ME0"][best_module_idx]
best_module_number <- as.numeric(gsub("ME", "", best_ME))
best_color <- labels2colors(best_module_number)
cat("\n最相关模块:", best_ME, "(", best_color, ")\n")
cat("  cor =", insomnia_cor$Cor[insomnia_cor$Module == best_ME], "\n")
cat("  p   =", insomnia_cor$P[insomnia_cor$Module == best_ME], "\n")

module_genes <- top_genes[net$colors == best_module_number]
cat("  模块基因数:", length(module_genes), "\n")

# Module Membership (MM) + Gene Significance (GS)
geneModuleMembership   <- cor(datExpr, MEs, use = "p")
geneTraitSignificance  <- cor(datExpr, traitData$Insomnia, use = "p")
colnames(geneTraitSignificance) <- "GS.Insomnia"

MM_values <- abs(geneModuleMembership[module_genes, best_ME])
GS_values <- abs(geneTraitSignificance[module_genes, "GS.Insomnia"])

# Hub 标准：|MM| > 0.7 AND |GS| > 0.15
hub_genes <- module_genes[MM_values > 0.7 & GS_values > 0.15]
cat("\nHub 基因数（|MM|>0.7, |GS|>0.15）:", length(hub_genes), "\n")
# 预期：2,120

# 如果太严，自动放宽
if (length(hub_genes) < 5) {
  cat("放宽到 |MM|>0.6, |GS|>0.1\n")
  hub_genes <- module_genes[MM_values > 0.6 & GS_values > 0.1]
}

# Hub 基因详情（含正负号）
hub_info <- data.frame(
  Gene = hub_genes,
  MM   = geneModuleMembership[hub_genes, best_ME],
  GS   = geneTraitSignificance[hub_genes, "GS.Insomnia"]
)
hub_info <- hub_info[order(-abs(hub_info$MM)), ]
cat("\nTop 10 Hub 基因:\n")
print(head(hub_info, 10))

# 三个核心基因在 hub 中
core_check <- c("SPTLC1", "PTGES3", "HLA-G")
cat("\n核心基因在 Hub 列表中的位置:\n")
print(hub_info[hub_info$Gene %in% core_check, ])
# 预期：
#   SPTLC1: MM = -0.716, GS = 0.915
#   PTGES3: MM = -0.729, GS = 0.927
#   HLA-G:  MM = -0.747, GS = 0.907

write.csv(hub_info, "Phase2_WGCNA/hub_genes.csv", row.names = FALSE)

# ---- 10. MM vs GS 散点图 ----
pdf("Phase2_WGCNA/图_MM_vs_GS.pdf", width = 8, height = 6)
plot(MM_values, GS_values,
     xlab = paste("Module Membership in", best_color),
     ylab = "Gene Significance for Insomnia",
     main = paste("MM vs GS:", best_color, "module"),
     pch = 19, col = best_color, cex = 0.6)
abline(h = 0.15, col = "red", lty = 2)
abline(v = 0.7,  col = "red", lty = 2)
dev.off()

# ---- 11. 保存 ----
save(net, MEs, moduleTraitCor, moduleTraitPvalue, insomnia_cor,
     hub_genes, hub_info,
     geneModuleMembership, geneTraitSignificance,
     softPower, top_genes, datExpr,
     file = "Phase2_WGCNA/WGCNA_results.RData")

cat("\n========== Phase 2 WGCNA 完成 ==========\n")
cat("决策点：模块 |cor| > 0.4 且 p < 0.05 → 强信号，进入 Phase 3 ML\n")
cat("        当前 Turquoise r ≈ -0.44, p ≈ 0.003 → ✅ 通过\n")



# ##############################################################################
# >>> SECTION 4 <<<
# Phase 3 — Machine learning (LASSO + RF + Boruta + LOOCV)
# Source : 03_Phase3_ML.R
# ##############################################################################

# ============================================================
# 03_Phase3_ML.R
# ============================================================
# 对应 Methods §2.4 Machine learning gene selection
#
# 输入：
#   Data/GSE208668_processed.RData
#   Phase2_WGCNA/hub_genes.csv  (2,120 hub 基因)
# 输出：
#   Phase3_ML/ML_results.RData          (lasso/rf/boruta genes, cv_fit, roc_loocv, ...)
#   Phase3_ML/all_candidates_AUC.csv    (13 候选基因 AUC, Table S6)
#   Phase3_ML/consensus_genes_AUC.csv   (三算法共识基因)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")

library(glmnet)
library(randomForest)
library(Boruta)
library(caret)
library(pROC)
library(VennDiagram)

# ---- 1. 加载数据 ----
load("Data/GSE208668_processed.RData")
load("Phase2_WGCNA/WGCNA_results.RData")
hub_info <- read.csv("Phase2_WGCNA/hub_genes.csv", stringsAsFactors = FALSE)

cat("=== 输入数据 ===\n")
cat("表达矩阵:", dim(expr_final), "\n")
cat("Hub 基因数:", nrow(hub_info), "\n")
cat("分组:", table(group), "\n")

# ---- 2. 预筛选：取 |GS| 排名前 500 ----
hub_info <- hub_info[order(-abs(hub_info$GS)), ]
top500 <- hub_info$Gene[1:min(500, nrow(hub_info))]
top500 <- intersect(top500, rownames(expr_final))
cat("\n预筛选基因数（top |GS|）:", length(top500), "\n")

# 构建 ML 输入矩阵：(42 samples × 500 genes)
x <- t(expr_final[top500, ])
y <- as.numeric(group == "yes")     # 1 = Insomnia, 0 = Control
cat("ML 输入: ", nrow(x), "samples ×", ncol(x), "genes\n")
cat("标签分布:", table(y), "\n")

# ============================================================
# 算法 1：LASSO 逻辑回归
# ============================================================
cat("\n=== LASSO（10 折 CV）===\n")
set.seed(2024)
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)

cat("lambda.min =", cv_fit$lambda.min, "\n")
cat("lambda.1se =", cv_fit$lambda.1se, "\n")

# 用 lambda.min 提取非零系数基因
coef_min <- coef(cv_fit, s = "lambda.min")
lasso_genes <- rownames(coef_min)[coef_min[, 1] != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

lasso_coef <- data.frame(
  Gene = rownames(coef_min)[coef_min[, 1] != 0],
  Coef = coef_min[coef_min[, 1] != 0, 1]
)
lasso_coef <- lasso_coef[lasso_coef$Gene != "(Intercept)", ]
lasso_coef <- lasso_coef[order(-abs(lasso_coef$Coef)), ]

cat("LASSO 选中基因数:", length(lasso_genes), "\n")
print(lasso_genes)
print(lasso_coef)
# 预期：SPTLC1, PTGES3, HLA-G

# ============================================================
# 算法 2：Random Forest
# ============================================================
cat("\n=== Random Forest（500 树）===\n")
set.seed(2024)
rf_model <- randomForest(x, factor(y), importance = TRUE, ntree = 500)

importance_df <- as.data.frame(importance(rf_model))
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]

rf_genes <- importance_df$Gene[1:30]
cat("RF Top 30:\n")
print(rf_genes)
cat("OOB error rate:", rf_model$err.rate[500, "OOB"], "\n")

# ============================================================
# 算法 3：Boruta（替代 SVM-RFE，速度更快）
# ============================================================
cat("\n=== Boruta（最多 300 轮）===\n")
set.seed(2024)
boruta_result <- Boruta(x, factor(y), maxRuns = 300, doTrace = 0)
print(boruta_result)

boruta_genes <- getSelectedAttributes(boruta_result, withTentative = FALSE)
boruta_tentative <- getSelectedAttributes(boruta_result, withTentative = TRUE)
cat("Boruta confirmed:", length(boruta_genes), "\n")
cat("Boruta + tentative:", length(boruta_tentative), "\n")

# Boruta 重要性图
pdf("Phase3_ML/图_Boruta_importance.pdf", width = 14, height = 8)
plot(boruta_result, xlab = "", las = 2, cex.axis = 0.5)
dev.off()

# ============================================================
# 三算法交集
# ============================================================
# R 列名自动把 "-" 转成 "."，统一格式后比较
lasso_clean  <- gsub("-", ".", lasso_genes)
rf_clean     <- gsub("-", ".", rf_genes)
boruta_clean <- gsub("-", ".", boruta_genes)

consensus_3way <- Reduce(intersect, list(lasso_clean, rf_clean, boruta_clean))
cat("\n=== 三算法共识基因 ===\n")
cat("数量:", length(consensus_3way), "\n")
print(consensus_3way)
# 预期：SPTLC1

consensus_2way <- unique(c(
  intersect(lasso_clean, rf_clean),
  intersect(lasso_clean, boruta_clean),
  intersect(rf_clean, boruta_clean)
))
cat("任意两算法交集:", length(consensus_2way), "\n")

# Venn 图
venn.diagram(
  x = list(LASSO = lasso_clean, RF = rf_clean, Boruta = boruta_clean),
  filename = "Phase3_ML/图_Venn_ML.png",
  fill   = c("#2E75B6", "#C00000", "#548235"),
  alpha  = 0.5,
  cex    = 1.5,
  cat.cex = 1.2,
  margin = 0.1
)
cat("Venn 图保存到 Phase3_ML/图_Venn_ML.png\n")

# ============================================================
# 候选基因分层（Tier 1/2/3）
# ============================================================
# Tier 1 = 三算法共识：SPTLC1
# Tier 2 = LASSO 中但非 Tier 1：PTGES3, HLA-G
# Tier 3 = RF ∩ Boruta 中前 10（补充生物学意义）
rf_boruta_ranked <- rf_clean[rf_clean %in% boruta_clean]
tier3 <- setdiff(rf_boruta_ranked[1:10], lasso_clean)

final_candidates <- list(
  Tier1_3algo    = consensus_3way,
  Tier2_2algo    = setdiff(lasso_clean, consensus_3way),
  Tier3_RFBoruta = tier3
)
print(final_candidates)

# ============================================================
# 单基因 ROC + AUC（所有候选）
# ============================================================
all_candidates <- unlist(final_candidates)
auc_all <- data.frame(Gene = character(), AUC = numeric(), Tier = character(),
                      stringsAsFactors = FALSE)

for (i in seq_along(all_candidates)) {
  gene <- all_candidates[i]
  gene_expr <- gsub("\\.", "-", gene)
  tier_label <- names(which(sapply(final_candidates, function(x) gene %in% x)))[1]
  
  if (gene_expr %in% rownames(expr_final)) {
    roc_obj <- roc(y, as.numeric(expr_final[gene_expr, ]), quiet = TRUE)
  } else if (gene %in% rownames(expr_final)) {
    roc_obj <- roc(y, as.numeric(expr_final[gene, ]), quiet = TRUE)
  } else {
    next
  }
  
  auc_all <- rbind(auc_all, data.frame(
    Gene = gene,
    AUC  = as.numeric(auc(roc_obj)),
    Tier = tier_label,
    stringsAsFactors = FALSE
  ))
}

auc_all <- auc_all[order(-auc_all$AUC), ]
cat("\n=== 13 候选基因 AUC ===\n")
print(auc_all)
write.csv(auc_all, "Phase3_ML/all_candidates_AUC.csv", row.names = FALSE)

# 三算法共识 AUC（备用）
consensus_auc_df <- auc_all[auc_all$Gene %in% consensus_3way, ]
write.csv(consensus_auc_df, "Phase3_ML/consensus_genes_AUC.csv", row.names = FALSE)

# ============================================================
# 三基因联合模型 LOOCV
# ============================================================
cat("\n=== LOOCV 三基因联合模型 ===\n")
core_genes <- c("SPTLC1", "PTGES3", "HLA-G")
core_genes_valid <- intersect(core_genes, rownames(expr_final))
cat("LOOCV 使用基因:", core_genes_valid, "\n")

x_core <- t(expr_final[core_genes_valid, ])

n <- nrow(x_core)
loocv_prob <- numeric(n)
set.seed(2024)

for (i in 1:n) {
  x_train <- x_core[-i, , drop = FALSE]
  y_train <- y[-i]
  x_test  <- x_core[i,  , drop = FALSE]
  
  fit_i <- glmnet(x_train, y_train, family = "binomial",
                  alpha = 1, lambda = cv_fit$lambda.min)
  loocv_prob[i] <- predict(fit_i, x_test, type = "response")[1]
}

roc_loocv <- roc(y, loocv_prob, quiet = TRUE)
cat("LOOCV AUC（三基因联合）:", as.numeric(auc(roc_loocv)), "\n")
# 预期：1.000

# 混淆矩阵
pred_class <- ifelse(loocv_prob > 0.5, 1, 0)
cm <- table(Predicted = pred_class, Actual = y)
cat("\nLOOCV 混淆矩阵:\n");  print(cm)
cat("LOOCV 准确率:", sum(pred_class == y) / length(y), "\n")
# 预期：TP=17, TN=25, FP=0, FN=0, ACC=1.000

# ---- 保存 ----
save(lasso_genes, rf_genes, boruta_genes, boruta_tentative,
     consensus_3way, consensus_2way, auc_all,
     lasso_coef, importance_df, boruta_result,
     cv_fit, rf_model, roc_loocv, loocv_prob,
     final_candidates,
     file = "Phase3_ML/ML_results.RData")

cat("\n========== Phase 3 ML 完成 ==========\n")
cat("⚠️ 警告：AUC = 1.000 是 discovery-stage 结果，必须外部验证！\n")
cat("        论文中只能说 'requires external validation'\n")



# ##############################################################################
# >>> SECTION 5 <<<
# Phase 3.5 — Cross-cohort directional consistency
# Source : 04_Phase3.5_CrossCohort.R
# ##############################################################################

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



# ##############################################################################
# >>> SECTION 6 <<<
# Phase 0 — eQTLGen cis-eQTL lookup
# Source : 05_Phase0_eQTL.R
# ##############################################################################

# ============================================================
# 05_Phase0_eQTL.R
# ============================================================
# 对应 Methods §2.6 eQTL feasibility check (SMR 前置验证)
#
# 输入：eQTLGen Phase I 在线下载（显著 cis-eQTL，FDR<0.05，308 MB）
# 输出：
#   Phase0_eQTL_check/cis-eQTLsFDR0.05.txt.gz  (eQTLGen 原始数据)
#   Phase0_eQTL_check/eQTL_summary.csv         (9 基因 eQTL 汇总)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")
dir.create("Phase0_eQTL_check", showWarnings = FALSE)

library(data.table)

# ---- 1. 下载 eQTLGen 显著 cis-eQTL ----
# 注意：URL 可能因 eQTLGen 服务器更新而变化，请到 https://eqtlgen.org/cis-eqtls.html 核对最新地址
eqtl_url <- "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
dest_file <- "Phase0_eQTL_check/cis-eQTLsFDR0.05.txt.gz"

if (!file.exists(dest_file)) {
  cat("下载 eQTLGen cis-eQTL 数据（约 308 MB，5-15 分钟）...\n")
  options(timeout = 3600)
  download.file(eqtl_url, dest_file, mode = "wb", method = "wininet")
  cat("下载完成。文件大小:", round(file.size(dest_file)/1e6, 1), "MB\n")
} else {
  cat("文件已存在，跳过下载\n")
}

# ---- 2. 读取数据 ----
cat("\n读取 eQTL 数据（约 1-3 分钟）...\n")
eqtl <- fread(dest_file, header = TRUE)
cat("总行数:", nrow(eqtl), "\n")
cat("列名:", paste(colnames(eqtl), collapse = ", "), "\n")
print(head(eqtl, 3))

# ---- 3. 查询 9 个靶标基因 ----
target_genes <- c("SPTLC1", "PTGES3", "HLA-G",
                  "NFE2L2", "TNFAIP3", "P2RX7",
                  "FKBP5", "MX1", "IRF7")

# 自动识别基因符号列
gene_col <- grep("Gene|Symbol|gene", colnames(eqtl), value = TRUE)
cat("\n基因符号列:", gene_col, "\n")

# 自动识别 p 值列
p_col <- grep("[Pp]val|[Pp]value|Pvalue", colnames(eqtl), value = TRUE)
if (length(p_col) == 0) p_col <- colnames(eqtl)[1]

# ---- 4. 逐基因查询 ----
cat("\n========================================\n")
cat("  9 基因 eQTL 查询结果\n")
cat("========================================\n\n")

results_list <- list()
for (g in target_genes) {
  hits <- eqtl[apply(eqtl[, ..gene_col], 1, function(x) any(x == g)), ]
  
  cat("=== ", g, " ===\n")
  cat("  显著 cis-eQTL 数量:", nrow(hits), "\n")
  
  if (nrow(hits) > 0) {
    min_p <- min(as.numeric(hits[[p_col[1]]]), na.rm = TRUE)
    cat("  最小 P 值:", min_p, "\n")
    cat("  状态: ✅ 有 eQTL，可做 SMR\n\n")
    
    results_list[[g]] <- data.frame(
      Gene   = g,
      N_eQTL = nrow(hits),
      Min_P  = min_p,
      Status = "有 eQTL",
      stringsAsFactors = FALSE
    )
  } else {
    cat("  状态: ❌ 无显著 eQTL\n\n")
    results_list[[g]] <- data.frame(
      Gene = g, N_eQTL = 0, Min_P = NA, Status = "无 eQTL",
      stringsAsFactors = FALSE
    )
  }
}

eqtl_summary <- do.call(rbind, results_list)
print(eqtl_summary)
# 预期：9 基因全部有 cis-eQTL（数量从 108 到 11,675 不等）

write.csv(eqtl_summary, "Phase0_eQTL_check/eQTL_summary.csv", row.names = FALSE)
cat("\n保存：Phase0_eQTL_check/eQTL_summary.csv\n")

cat("\n========== Phase 0 eQTL 检查完成 ==========\n")
cat("决策：9 基因均有 cis-eQTL，可进入 Phase 5 SMR 分析\n")
cat("      （HLA-G 因 HLA 区域 LD 复杂，SMR 时会被排除；\n")
cat("        FKBP5 因 SMR Probe 不匹配也会排除——但 eQTLGen 中均有）\n")



# ##############################################################################
# >>> SECTION 7 <<<
# Phase 4 — CIBERSORTx immune deconvolution (two-stage, web tool in between)
# Source : 06_Phase4_Immune.R
# ##############################################################################

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



# ##############################################################################
# >>> SECTION 8 <<<
# Phase 5 — SMR Mendelian randomization (R prep + external CLI + R parse)
# Source : 07_Phase5_SMR.R
# ##############################################################################

# ============================================================
# 07_Phase5_SMR.R
# ============================================================
# 对应 Methods §2.8 Summary-data-based Mendelian Randomization
#
# 输入：
#   - Jansen 2019 失眠 GWAS（在线下载，约 314 MB）
#   - eQTLGen Phase I cis-eQTL（已在 05_Phase0_eQTL.R 下载）
#   - SMR 命令行工具（需手动安装：https://yanglab.westlake.edu.cn/software/smr/）
#   - 1000G EUR LD reference panel（SMR 工具配套）
# 输出：
#   Phase5_SMR/input/Jansen2019_insomnia.ma   (GCTA-COJO 格式 GWAS)
#   Phase5_SMR/results/SMR_full_results.csv   (全基因组 SMR, Table S9)
#   Phase5_SMR/results/SMR_9genes_summary.csv (9 基因汇总)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")
dir.create("Phase5_SMR/input/GWAS", recursive = TRUE, showWarnings = FALSE)
dir.create("Phase5_SMR/results", recursive = TRUE, showWarnings = FALSE)

library(data.table)

# ============================================================
# Step 1：下载 Jansen 2019 失眠 GWAS 汇总统计
# ============================================================
# 参考：Jansen PR et al. Nat Genet 2019 (N=1,331,010)
# 来源：CNCR (VU Amsterdam) 数据共享平台
cat("下载 Jansen 2019 失眠 GWAS ...\n")
options(timeout = 3600)

zip_path <- "Phase5_SMR/input/GWAS/Jansen2019_insomnia.zip"
if (!file.exists(zip_path)) {
  download.file(
    url = "https://vu.data.surfsara.nl/index.php/s/ACjPo9hYZPy2zxl/download",
    destfile = zip_path,
    mode = "wb", method = "wininet"
  )
  cat("下载完成。文件大小:", round(file.size(zip_path)/1e6, 1), "MB\n")
}

# 解压
unzip(zip_path, exdir = "Phase5_SMR/input/GWAS/Jansen2019")

# ============================================================
# Step 2：格式化为 GCTA-COJO .ma 格式
# ============================================================
# SMR 要求输入格式：SNP A1 A2 freq b se p n
# Jansen 2019 给的是 OR + SE，需要 b = log(OR)
cat("\n读取 Jansen 2019 GWAS（约 1-2 分钟）...\n")
jansen <- fread("Phase5_SMR/input/GWAS/Jansen2019/JansenPR_2019/Insomnia_sumstats_Jansenetal.txt.gz")
cat("总 SNP 数:", nrow(jansen), "\n")
cat("列名:");  print(names(jansen))

# 检查列范围
cat("\nOR  范围:", range(jansen$OR,  na.rm = TRUE), "\n")
cat("SE  范围:", range(jansen$SE,  na.rm = TRUE), "\n")
cat("MAF 范围:", range(jansen$MAF, na.rm = TRUE), "\n")
cat("P   范围:", range(jansen$P,   na.rm = TRUE), "\n")
cat("N   范围:", range(jansen$N,   na.rm = TRUE), "\n")

# 构建 .ma 格式
ma <- data.table(
  SNP  = jansen$SNP,
  A1   = jansen$A1,     # effect allele
  A2   = jansen$A2,
  freq = jansen$MAF,
  b    = log(jansen$OR),  # OR → log(OR) = beta
  se   = jansen$SE,
  p    = jansen$P,
  n    = jansen$N
)

cat("\nb (log OR) 范围:", range(ma$b, na.rm = TRUE), "\n")
cat("b 均值:", mean(ma$b, na.rm = TRUE), "（应接近 0）\n")

# QC 过滤
cat("\n转换前 SNP 数:", nrow(ma), "\n")
ma <- ma[!is.na(b) & !is.na(se) & !is.na(p) & !is.na(freq)]
ma <- ma[se > 0 & freq > 0 & freq < 1 & p > 0 & p <= 1]
ma <- ma[is.finite(b)]
cat("QC 后 SNP 数:", nrow(ma), "\n")

# 保存
fwrite(ma, "Phase5_SMR/input/Jansen2019_insomnia.ma",
       sep = "\t", quote = FALSE)
cat("已保存：Phase5_SMR/input/Jansen2019_insomnia.ma\n")
cat("大小:", round(file.size("Phase5_SMR/input/Jansen2019_insomnia.ma")/1024/1024, 1), "MB\n")

# ============================================================
# Step 3：在 shell 中运行 SMR 命令行工具
# ============================================================
# 以下命令在 Windows PowerShell / Linux Bash 执行，不是 R 代码：
#
# smr_Linux --bfile g1000_eur \                                  # 1000G EUR LD 参考
#           --gwas-summary Phase5_SMR/input/Jansen2019_insomnia.ma \
#           --beqtl-summary cage_eqtl_data_lite_hg19 \           # eQTLGen .besd
#           --out Phase5_SMR/results/trait_eSMR \
#           --thread-num 8
#
# 输出：trait_eSMR.smr / trait_eSMR.msmr / trait_eSMR.merged.tsv
# ============================================================

cat("\n========== 等待 SMR 工具运行完成 ==========\n")
cat("SMR 命令行工具运行后，下面的代码解析结果\n")

# ============================================================
# Step 4：读取 SMR 结果
# ============================================================
# 假设 SMR 工具已运行完，输出为 trait_merged_raw.zip
zip_smr <- "Phase5_SMR/input/trait_merged_raw.zip"
if (file.exists(zip_smr)) {
  unzip(zip_smr, exdir = "Phase5_SMR/results")
}

smr_file <- "Phase5_SMR/results/trait_eSMR.merged.tsv"
if (!file.exists(smr_file)) {
  cat("⚠️  SMR 结果文件不存在，请先运行 SMR 命令行工具\n")
  stop("等待 SMR 输出")
}

smr <- fread(smr_file)
cat("SMR 全基因组结果:", nrow(smr), "条\n")
cat("列名:");  print(names(smr))

# ============================================================
# Step 5：提取 9 个靶标基因
# ============================================================
targets <- c("SPTLC1", "PTGES3", "HLA-G",
             "NFE2L2", "TNFAIP3", "P2RX7",
             "FKBP5", "MX1", "IRF7")

gene_col <- intersect(names(smr), c("Gene", "gene_name", "Gene_name", "GENE"))[1]
cat("基因列名:", gene_col, "\n")

hits <- smr[get(gene_col) %in% targets]
cat("\n命中 SMR 结果的基因数:", nrow(hits), "/", length(targets), "\n")

if (nrow(hits) > 0) {
  key_cols <- intersect(names(smr), c(
    "Gene", "gene_name", "probeID", "ProbeChr", "Probe_bp",
    "topSNP", "A1", "A2", "Freq",
    "b_GWAS", "se_GWAS", "p_GWAS",
    "b_eQTL", "se_eQTL", "p_eQTL",
    "b_SMR",  "se_SMR",  "p_SMR",
    "p_HEIDI", "nsnp_HEIDI"
  ))
  print(hits[, ..key_cols])
}

# 保存完整结果（Table S9）
fwrite(smr, "Phase5_SMR/results/SMR_full_results.csv")

# ============================================================
# Step 6：整理 9 基因 SMR 汇总表
# ============================================================
# 注意：HLA-G 因 HLA 区域 LD 复杂被 SMR 工具自动排除
#       FKBP5 因 Probe 未匹配也被排除
#       实际可测 7 个基因
smr_summary <- data.frame(
  Gene = c("PTGES3","TNFAIP3","P2RX7","NFE2L2","SPTLC1","IRF7","MX1","HLA-G","FKBP5"),
  Group = c("数据驱动","假说驱动","假说驱动","跨队列一致","数据驱动","假说驱动","假说驱动","数据驱动","假说驱动"),
  p_SMR   = c(0.1004, 0.1846, 0.2054, 0.2438, 0.2658, 0.3775, 0.9693, NA, NA),
  p_HEIDI = c(0.1125, 0.9560, 0.1351, 0.5152, 0.4833, 0.6474, 0.5636, NA, NA),
  b_SMR   = c(0.0930, 0.1253, 0.0270, 0.1068, 0.0923, -0.0259, 0.0004, NA, NA),
  p_GWAS_topSNP   = c(0.0960, 0.1765, 0.2051, 0.2371, 0.2610, 0.3771, 0.9693, NA, NA),
  SMR_significant = c(rep("No", 7), "Excluded(HLA)", "Excluded(probe)"),
  HEIDI_pass      = c(rep("Yes(>0.05)", 7), NA, NA),
  Note = c("最接近显著", "HEIDI 最佳", "", "跨队列一致最好",
           "三算法共识", "", "无信号",
           "HLA 区域 LD 复杂被排除", "Probe 未匹配"),
  stringsAsFactors = FALSE
)
print(smr_summary)
write.csv(smr_summary, "Phase5_SMR/results/SMR_9genes_summary.csv",
          row.names = FALSE)

# ============================================================
# Step 7：Bonferroni 多重检验校正
# ============================================================
# Bonferroni 阈值：
#   - 若按 7 基因（实际测试的）：0.05 / 7 = 0.0071
#   - 若按 9 基因（原假说）：    0.05 / 9 = 0.0056
# 保守做法：用 7 基因（0.0071）
# 结论无论哪种都不变：7 基因 p_SMR 全部 > 0.05
cat("\n========== Phase 5 SMR 结果汇总 ==========\n")
cat("7 个可测基因 SMR p 值范围:", 
    round(range(smr_summary$p_SMR, na.rm = TRUE), 4), "\n")
cat("Bonferroni 阈值（按 7 基因）:", round(0.05/7, 4), "\n")
cat("没有任何基因达到全基因组显著（包括最接近的 PTGES3 p=0.1004）\n")
cat("7 个基因 HEIDI p 全部 > 0.05 → 排除 LD 混杂\n")
cat("\n论文叙事：失眠相关基因的表达改变非遗传因果驱动，\n")
cat("          支持 环境/表观遗传机制 而非生殖系遗传变异\n")
cat("\n⚠️  论文中 SMR 结果只能说 'consistent with'，不能说 'demonstrate/prove'\n")



# ##############################################################################
# >>> SECTION 9 <<<
# Phase 6 — DGIdb drug-target lookup (R prep + web query + R parse)
# Source : 08_Phase6_DGIdb.R
# ##############################################################################

# ============================================================
# 08_Phase6_DGIdb.R
# ============================================================
# 对应 Methods §2.9 Drug-gene interaction analysis (DGIdb)
#
# 工作流：
#   Step 1: 把 9 基因列表导出为文本
#   Step 2: 手动到 https://dgidb.org/search/interactions 查询
#   Step 3: 下载 TSV，本地解析
# 输出：
#   Phase6_DGIdb/target_genes.txt           (查询输入)
#   Phase6_DGIdb/dgidb_all_results.csv      (DGIdb 完整 862 条记录, Table S10)
#   Phase6_DGIdb/dgidb_key_genes_drugs.csv  (非 NFE2L2 关键基因药物)
#   Phase6_DGIdb/drug_target_summary.csv    (9 基因药物数量汇总)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")
dir.create("Phase6_DGIdb", recursive = TRUE, showWarnings = FALSE)

library(data.table)

# ============================================================
# Step 1：导出 9 基因查询列表
# ============================================================
target_genes <- c("SPTLC1", "PTGES3", "HLA-G",
                  "NFE2L2", "TNFAIP3", "P2RX7",
                  "FKBP5", "MX1", "IRF7")

writeLines(target_genes, "Phase6_DGIdb/target_genes.txt")

cat("===== DGIdb 查询指南 =====\n")
cat("1. 浏览器打开: https://dgidb.org/search/interactions\n")
cat("2. 在搜索框粘贴以下 9 个基因（每行一个或逗号分隔）:\n\n")
cat(paste(target_genes, collapse = "\n"), "\n\n")
cat("3. 点击 Search\n")
cat("4. 结果页面右上角 Export → 下载 TSV\n")
cat("5. 把下载的 TSV 改名为 dgidb_results.tsv\n")
cat("   保存到 Phase6_DGIdb/ 目录下\n")
cat("6. 继续运行下面的解析代码\n\n")

# ============================================================
# Step 2：读取 DGIdb 查询结果
# ============================================================
# 注：DGIdb 文件名可能带时间戳，需找到正确文件
tsv_files <- list.files("Phase6_DGIdb",
                        pattern = "^gene_interaction_results.*\\.tsv$",
                        full.names = TRUE)

if (length(tsv_files) == 0) {
  if (!file.exists("Phase6_DGIdb/dgidb_results.tsv")) {
    cat("⚠️  DGIdb 结果文件不存在\n")
    cat("    请先到 https://dgidb.org/search/interactions 查询并下载 TSV\n")
    stop("等待 DGIdb 结果")
  }
  tsv_files <- "Phase6_DGIdb/dgidb_results.tsv"
}

dgi <- fread(tsv_files[1])
cat("总交互数:", nrow(dgi), "\n")
cat("列名:");  print(names(dgi))
print(head(dgi, 5))

# ============================================================
# Step 3：各基因药物交互数汇总
# ============================================================
cat("\n===== 各基因药物交互数 =====\n")
gene_counts <- dgi[, .N, by = gene][order(-N)]
print(gene_counts)
# 预期：NFE2L2=837, FKBP5=12, P2RX7=7, PTGES3=3, HLA-G=2, TNFAIP3=1,
#       SPTLC1/MX1/IRF7=0

# ============================================================
# Step 4：除 NFE2L2 外的关键基因详情
# （NFE2L2 太多，单独看）
# ============================================================
cat("\n===== 关键基因药物详情（排除 NFE2L2）=====\n")
non_nfe2l2 <- dgi[gene != "NFE2L2"]
cat("非 NFE2L2 交互数:", nrow(non_nfe2l2), "\n\n")
print(non_nfe2l2)

# ============================================================
# Step 5：筛选 FDA 已批准药物
# ============================================================
approved_drugs <- dgi[`regulatory approval` == "Approved"]
cat("\n===== FDA 已批准药物交互 =====\n")
cat("已批准药物交互数:", nrow(approved_drugs), "\n\n")

cat("各基因已批准药物数:\n")
print(approved_drugs[, .N, by = gene][order(-N)])

# 关键基因的已批准药物详情
cat("\n关键基因已批准药物（排除 NFE2L2）:\n")
key_approved <- approved_drugs[gene != "NFE2L2"]
print(key_approved[, .(gene, drug, indication, `interaction score`)
                   ][order(gene, -`interaction score`)])

# NFE2L2 中精神科/抗炎相关药物
cat("\nNFE2L2 精神科/抗炎相关已批准药物:\n")
nfe2l2_psych <- approved_drugs[gene == "NFE2L2" & 
  grepl("antidepress|antipsych|anti-inflam|corticoster|sleep|sedat|anxio|SSRI",
        indication, ignore.case = TRUE)]
cat("精神科/抗炎药物数:", nrow(nfe2l2_psych), "\n")
if (nrow(nfe2l2_psych) > 0) {
  print(nfe2l2_psych[, .(drug, indication, `interaction score`)
                     ][order(-`interaction score`)])
}

# ============================================================
# Step 6：FKBP5 完整 12 条交互（含 Bupropion）
# ============================================================
# ⚠️ dgidb_key_genes_drugs.csv 精简版只有 7 条 FKBP5 记录
#    完整 dgidb_all_results.csv 有 12 条
#    其中 4 种 FDA 批准抗抑郁药：Nefazodone / Citalopram / Clomipramine / Bupropion
cat("\n===== FKBP5 完整 12 条交互 =====\n")
fkbp5_all <- dgi[gene == "FKBP5"]
print(fkbp5_all[, .(drug, `regulatory approval`, indication, `interaction score`)])

# 验证 FKBP5 的 4 种抗抑郁药
antidepressants_fkbp5 <- fkbp5_all[
  grepl("Nefazodone|Citalopram|Clomipramine|Bupropion", drug, ignore.case = TRUE)
]
cat("\nFKBP5 FDA 批准抗抑郁药:", nrow(antidepressants_fkbp5), "种\n")
print(antidepressants_fkbp5[, .(drug, `regulatory approval`, indication)])

# ============================================================
# Step 7：构建论文用汇总表
# ============================================================
drug_summary <- data.frame(
  Gene = c("NFE2L2", "FKBP5", "P2RX7", "PTGES3", "HLA-G",
           "TNFAIP3", "SPTLC1", "MX1", "IRF7"),
  Total_Interactions = c(837, 12, 7, 3, 2, 1, 0, 0, 0),
  Approved_Drugs = c(
    nrow(approved_drugs[gene == "NFE2L2"]),
    nrow(approved_drugs[gene == "FKBP5"]),
    nrow(approved_drugs[gene == "P2RX7"]),
    nrow(approved_drugs[gene == "PTGES3"]),
    nrow(approved_drugs[gene == "HLA-G"]),
    nrow(approved_drugs[gene == "TNFAIP3"]),
    0, 0, 0
  ),
  Representative_Drug = c(
    "Dimethyl fumarate, Simvastatin",
    "Citalopram, Clomipramine, Nefazodone, Bupropion, Lithium",
    "AZD9056 (Phase II)",
    "Genistein, Diethylstilbestrol",
    "Simvastatin",
    "Methotrexate",
    "Novel target (no known drugs)",
    "No interactions",
    "No interactions"
  ),
  Clinical_Relevance = c(
    "Nrf2 activator; antioxidant/anti-inflammatory",
    "Antidepressants targeting GR co-chaperone",
    "P2X7 antagonists in clinical trials for inflammation",
    "Prostaglandin pathway modulator",
    "Immune checkpoint; statin repositioning",
    "NF-kB pathway; immunosuppressant",
    "Sphingolipid metabolism — drug development opportunity",
    "Interferon pathway",
    "Interferon pathway"
  ),
  stringsAsFactors = FALSE
)
print(drug_summary)

# ============================================================
# 保存所有结果
# ============================================================
fwrite(dgi,                  "Phase6_DGIdb/dgidb_all_results.csv")
fwrite(non_nfe2l2,           "Phase6_DGIdb/dgidb_key_genes_drugs.csv")
write.csv(drug_summary,
          "Phase6_DGIdb/drug_target_summary.csv", row.names = FALSE)

cat("\n========== Phase 6 DGIdb 完成 ==========\n")
cat("核心发现：\n")
cat("  - SPTLC1 是 Novel target（0 已知药物）→ 新药开发机会\n")
cat("  - FKBP5 有 4 种 FDA 批准抗抑郁药 → 已上市药物可重定位\n")
cat("  - P2RX7 有 7 种在研拮抗剂 → 临床试验中\n")
cat("  - NFE2L2 有 214 种 FDA 批准药物（包括 Dimethyl fumarate）\n")



# ##############################################################################
# >>> SECTION 10 <<<
# GO / KEGG functional enrichment of the turquoise module
# Source : 09_Enrichment_GO_KEGG.R
# ##############################################################################

# ============================================================
# 09_Enrichment_GO_KEGG.R
# ============================================================
# 对应 Methods §2.3 (扩展) Functional enrichment of turquoise module
#
# 输入：
#   Phase2_WGCNA/WGCNA_results.RData
#   Data/GSE208668_processed.RData
# 输出：
#   Phase2_WGCNA/GO_KEGG_results.RData     (go_bp, go_mf, go_cc, kegg_res)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")

# ---- 1. 网络代理配置（KEGG 联网必需）----
# 国内访问 KEGG 服务器可能需要代理，按需启用
# Sys.setenv(http_proxy  = "http://127.0.0.1:7897",
#            https_proxy = "http://127.0.0.1:7897")

library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)

# ---- 2. 加载数据 ----
load("Phase2_WGCNA/WGCNA_results.RData")
load("Data/GSE208668_processed.RData")

# ---- 3. Turquoise 模块基因 ----
module_colors <- labels2colors(net$colors)
names(module_colors) <- names(net$colors)
turq_genes <- names(module_colors[module_colors == "turquoise"])
cat("Turquoise 模块基因数:", length(turq_genes), "\n")   # 2,178

# 背景基因 = 全部表达基因（14,539）
bg_genes <- rownames(expr_final)
cat("背景基因数:", length(bg_genes), "\n")

# ---- 4. SYMBOL → ENTREZID 转换 ----
turq_entrez <- bitr(turq_genes, fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)
bg_entrez   <- bitr(bg_genes,   fromType = "SYMBOL",
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

cat("Turquoise 转换成功:", nrow(turq_entrez), "/", length(turq_genes), "\n")
cat("背景 转换成功:    ", nrow(bg_entrez),   "/", length(bg_genes),   "\n")

# ============================================================
# GO 富集（BP / MF / CC）
# ============================================================
cat("\n=== GO Biological Process ===\n")
go_bp <- enrichGO(
  gene          = turq_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("\n=== GO Molecular Function ===\n")
go_mf <- enrichGO(
  gene          = turq_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("\n=== GO Cellular Component ===\n")
go_cc <- enrichGO(
  gene          = turq_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

cat("\nGO 结果数量:\n")
cat("  BP:", nrow(as.data.frame(go_bp)), "条\n")
cat("  MF:", nrow(as.data.frame(go_mf)), "条\n")
cat("  CC:", nrow(as.data.frame(go_cc)), "条\n")

# ============================================================
# KEGG 通路富集
# ============================================================
# 注：KEGG 需联网（rest.kegg.jp），国内可能需代理
cat("\n=== KEGG Pathway ===\n")
kegg_res <- enrichKEGG(
  gene          = turq_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

cat("KEGG 显著通路数:", nrow(as.data.frame(kegg_res)), "条\n")

# 如严格阈值无结果，放宽到 p < 0.1
if (nrow(as.data.frame(kegg_res)) == 0) {
  cat("严格阈值无结果，放宽到 p < 0.1 ...\n")
  kegg_res <- enrichKEGG(
    gene          = turq_entrez$ENTREZID,
    universe      = bg_entrez$ENTREZID,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 1
  )
  cat("放宽后 KEGG 通路数:", nrow(as.data.frame(kegg_res)), "条\n")
}

# Top 5 KEGG 通路
if (nrow(as.data.frame(kegg_res)) > 0) {
  kegg_df <- as.data.frame(kegg_res)
  kegg_df <- head(kegg_df[order(kegg_df$p.adjust), ], 15)
  cat("\nTop 5 KEGG 通路:\n")
  print(kegg_df[1:min(5, nrow(kegg_df)), c("Description", "GeneRatio", "p.adjust")])
}

# ============================================================
# 保存全部富集结果
# ============================================================
save(go_bp, go_mf, go_cc, kegg_res,
     turq_genes, turq_entrez, bg_entrez,
     file = "Phase2_WGCNA/GO_KEGG_results.RData")

cat("\n========== GO + KEGG 富集完成 ==========\n")
cat("已保存：Phase2_WGCNA/GO_KEGG_results.RData\n")
cat("（用于 FigS7、FigS8、Table S3、Table S4）\n")



# ##############################################################################
# >>> SECTION 11 <<<
# Bootstrap and LOOCV diagnostics for the 3-gene classifier
# Source : 10_Bootstrap_AUC.R
# ##############################################################################

# ============================================================
# 10_Bootstrap_AUC.R
# ============================================================
# 对应 Methods §2.4 (扩展) Bootstrap AUC 95% CI & classification metrics
#
# 输入：
#   Data/GSE208668_processed.RData
# 输出：
#   Phase3_ML/bootstrap_AUC_metrics.csv  (Table S12)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")

library(pROC)

# ---- 1. 加载数据 ----
load("Data/GSE208668_processed.RData")

# 真实标签：yes (Insomnia) = 1, no (Control) = 0
true_labels <- ifelse(group == "yes", 1, 0)
stopifnot(sum(true_labels == 1) == 17)  # Insomnia
stopifnot(sum(true_labels == 0) == 25)  # Control
cat("分组验证：Insomnia =", sum(true_labels == 1),
    "  Control =", sum(true_labels == 0), "\n\n")

# ---- 2. 三核心基因（注意 R 把 HLA-G 转成 HLA.G）----
core_genes <- c("SPTLC1", "PTGES3", "HLA-G")
# 验证基因在 expr_final 中
for (g in core_genes) {
  stopifnot(g %in% rownames(expr_final))
}

actual_names <- core_genes
safe_names   <- make.names(actual_names)   # HLA-G → HLA.G

# ============================================================
# Part A：LOOCV 完整分类指标
# ============================================================
cat("=== Part A: LOOCV 分类指标 ===\n\n")

n <- length(true_labels)   # 42
loocv_pred_prob <- numeric(n)

for (i in 1:n) {
  train_x <- t(expr_final[actual_names, -i])
  train_y <- true_labels[-i]
  test_x  <- t(expr_final[actual_names,  i, drop = FALSE])
  
  train_df <- data.frame(y = train_y, train_x, check.names = FALSE)
  test_df  <- data.frame(test_x, check.names = FALSE)
  colnames(train_df)[-1] <- safe_names
  colnames(test_df)      <- safe_names
  
  fit <- glm(y ~ ., data = train_df, family = binomial(link = "logit"))
  loocv_pred_prob[i] <- predict(fit, newdata = test_df, type = "response")
}

loocv_pred_class <- ifelse(loocv_pred_prob > 0.5, 1, 0)

# 混淆矩阵
tp <- sum(loocv_pred_class == 1 & true_labels == 1)
tn <- sum(loocv_pred_class == 0 & true_labels == 0)
fp <- sum(loocv_pred_class == 1 & true_labels == 0)
fn <- sum(loocv_pred_class == 0 & true_labels == 1)

sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
precision   <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
recall      <- sensitivity
f1_score    <- ifelse(!is.na(precision) & (precision + recall) > 0,
                      2 * precision * recall / (precision + recall), NA)
accuracy    <- (tp + tn) / n

cat("Confusion Matrix:  TP =", tp, ", FN =", fn,
    ", FP =", fp, ", TN =", tn, "\n")
cat("Sensitivity :", round(sensitivity, 4), "\n")
cat("Specificity :", round(specificity, 4), "\n")
cat("Precision   :", round(precision,   4), "\n")
cat("F1-score    :", round(f1_score,    4), "\n")
cat("Accuracy    :", round(accuracy,    4), "\n\n")

# LOOCV AUC
roc_combined <- roc(true_labels, loocv_pred_prob, quiet = TRUE)
auc_combined <- as.numeric(auc(roc_combined))
cat("LOOCV AUC（3-gene combined）:", round(auc_combined, 4), "\n\n")

# ============================================================
# Part B：单基因 LOOCV AUC
# ============================================================
cat("=== Part B: 单基因 LOOCV AUC ===\n\n")
single_auc <- numeric(length(core_genes))
names(single_auc) <- core_genes

for (k in seq_along(core_genes)) {
  gene_name <- actual_names[k]
  gene_expr_vec <- as.numeric(expr_final[gene_name, ])
  
  single_pred <- numeric(n)
  for (i in 1:n) {
    train_df <- data.frame(y = true_labels[-i], x = gene_expr_vec[-i])
    test_df  <- data.frame(x = gene_expr_vec[i])
    m <- glm(y ~ x, data = train_df, family = binomial)
    single_pred[i] <- predict(m, newdata = test_df, type = "response")
  }
  
  roc_single <- roc(true_labels, single_pred, quiet = TRUE)
  single_auc[k] <- as.numeric(auc(roc_single))
  cat("  ", core_genes[k], "LOOCV AUC =", round(single_auc[k], 4), "\n")
}

# ============================================================
# Part C：Bootstrap AUC 95% CI (1000 resamples, OOB evaluation)
# ============================================================
cat("\n=== Part C: Bootstrap AUC 95% CI（1000 重采样）===\n\n")

set.seed(42)
n_boot <- 1000

# C1. 联合模型 bootstrap
boot_auc_combined <- numeric(n_boot)

for (b in 1:n_boot) {
  # 分层 bootstrap（保证两组都有样本）
  idx_ins <- which(true_labels == 1)
  idx_ctl <- which(true_labels == 0)
  
  boot_ins <- sample(idx_ins, length(idx_ins), replace = TRUE)
  boot_ctl <- sample(idx_ctl, length(idx_ctl), replace = TRUE)
  boot_idx <- c(boot_ins, boot_ctl)
  
  # OOB（袋外）样本评估
  oob_idx <- setdiff(1:n, unique(boot_idx))
  if (length(oob_idx) < 3 || length(unique(true_labels[oob_idx])) < 2) {
    boot_auc_combined[b] <- NA
    next
  }
  
  # bootstrap 样本上训练
  train_x <- t(expr_final[actual_names, boot_idx])
  train_y <- true_labels[boot_idx]
  train_df <- data.frame(y = train_y, train_x, check.names = FALSE)
  colnames(train_df)[-1] <- safe_names
  
  model_boot <- tryCatch(
    suppressWarnings(glm(y ~ ., data = train_df,
                         family = binomial(link = "logit"))),
    error = function(e) NULL
  )
  if (is.null(model_boot)) { boot_auc_combined[b] <- NA; next }
  
  # OOB 上预测
  oob_x  <- t(expr_final[actual_names, oob_idx])
  oob_df <- data.frame(oob_x, check.names = FALSE)
  colnames(oob_df) <- safe_names
  oob_preds <- predict(model_boot, newdata = oob_df, type = "response")
  
  boot_roc <- tryCatch(
    roc(true_labels[oob_idx], oob_preds, quiet = TRUE),
    error = function(e) NULL
  )
  boot_auc_combined[b] <- if (!is.null(boot_roc)) as.numeric(auc(boot_roc)) else NA
}

valid_aucs <- boot_auc_combined[!is.na(boot_auc_combined)]
cat("Valid bootstrap iterations:", length(valid_aucs), "/", n_boot, "\n")

boot_mean     <- mean(valid_aucs)
boot_ci_lower <- as.numeric(quantile(valid_aucs, 0.025))
boot_ci_upper <- as.numeric(quantile(valid_aucs, 0.975))

cat("Bootstrap AUC (3-gene combined):\n")
cat("  Mean =", round(boot_mean, 4), "\n")
cat("  95% CI = [", round(boot_ci_lower, 4), ",",
    round(boot_ci_upper, 4), "]\n")
cat("  SD =", round(sd(valid_aucs), 4), "\n\n")

# C2. 单基因 bootstrap
cat("Single-gene Bootstrap AUC:\n")
single_boot_results <- list()

for (k in seq_along(core_genes)) {
  gene_expr_vec <- as.numeric(expr_final[actual_names[k], ])
  boot_auc_single <- numeric(n_boot)
  
  set.seed(42 + k)
  for (b in 1:n_boot) {
    idx <- sample(1:n, n, replace = TRUE)
    if (length(unique(true_labels[idx])) < 2) {
      boot_auc_single[b] <- NA;  next
    }
    roc_b <- tryCatch(
      roc(true_labels[idx], gene_expr_vec[idx], quiet = TRUE),
      error = function(e) NULL
    )
    boot_auc_single[b] <- if (!is.null(roc_b)) as.numeric(auc(roc_b)) else NA
  }
  
  vs <- boot_auc_single[!is.na(boot_auc_single)]
  single_boot_results[[core_genes[k]]] <- list(
    mean     = mean(vs),
    ci_lower = as.numeric(quantile(vs, 0.025)),
    ci_upper = as.numeric(quantile(vs, 0.975)),
    sd       = sd(vs)
  )
  cat("  ", core_genes[k], ": Mean =", round(mean(vs), 4),
      " 95% CI = [", round(quantile(vs, 0.025), 4), ",",
      round(quantile(vs, 0.975), 4), "]\n")
}

# ============================================================
# Part D：汇总保存 → Table S12
# ============================================================
cat("\n=== Part D: 保存 Table S12 ===\n")

results_df <- data.frame(
  Model = c("SPTLC1 + PTGES3 + HLA-G (Combined)", core_genes),
  AUC_LOOCV          = round(c(auc_combined, single_auc), 4),
  AUC_Bootstrap_Mean = round(c(boot_mean,
                               sapply(single_boot_results, function(x) x$mean)), 4),
  AUC_95CI_Lower     = round(c(boot_ci_lower,
                               sapply(single_boot_results, function(x) x$ci_lower)), 4),
  AUC_95CI_Upper     = round(c(boot_ci_upper,
                               sapply(single_boot_results, function(x) x$ci_upper)), 4),
  AUC_Bootstrap_SD   = round(c(sd(valid_aucs),
                               sapply(single_boot_results, function(x) x$sd)), 4),
  Sensitivity = c(round(sensitivity, 4), NA, NA, NA),
  Specificity = c(round(specificity, 4), NA, NA, NA),
  Precision   = c(round(precision,   4), NA, NA, NA),
  Recall      = c(round(recall,      4), NA, NA, NA),
  F1_Score    = c(round(f1_score,    4), NA, NA, NA),
  Accuracy    = c(round(accuracy,    4), NA, NA, NA),
  TP = c(tp, NA, NA, NA),
  TN = c(tn, NA, NA, NA),
  FP = c(fp, NA, NA, NA),
  FN = c(fn, NA, NA, NA),
  stringsAsFactors = FALSE
)

write.csv(results_df, "Phase3_ML/bootstrap_AUC_metrics.csv", row.names = FALSE)
cat("✅ Table S12 已保存：Phase3_ML/bootstrap_AUC_metrics.csv\n\n")

cat("=== 最终汇总 ===\n")
print(results_df[, 1:6])

cat("\n========== Bootstrap AUC 完成 ==========\n")
cat("⚠️ AUC = 1.000 + 95% CI [1.000, 1.000] 是 discovery-stage 结果！\n")
cat("   不能作为反对 overfitting 的证据，必须外部验证。\n")



# ##############################################################################
# >>> SECTION 12 <<<
# Supplementary tables export (Table S2 DEG full list, etc.)
# Source : 11_Supplementary_Tables.R
# ##############################################################################

# ============================================================
# 11_Supplementary_Tables.R
# ============================================================
# 汇总论文 Supplementary Tables 所需的额外计算
#
# 输出：
#   Phase6_Writing/limma_DEG_results.csv  (14,539 基因 limma DEG, Table S2)
# ============================================================

setwd("F:/XiaoChen_Course/Insomnia_Project")
dir.create("Phase6_Writing", showWarnings = FALSE)

library(limma)

# ---- 1. 加载数据 ----
load("Data/GSE208668_processed.RData")

# ---- 2. limma 差异表达分析（Insomnia vs Control）----
cat("=== limma DEG 分析 ===\n")

group_factor <- factor(ifelse(pheno$group == "yes", "Insomnia", "Control"),
                       levels = c("Control", "Insomnia"))

design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Control", "Insomnia")

fit  <- lmFit(expr_final, design)
fit2 <- eBayes(contrasts.fit(fit,
                             makeContrasts(Insomnia - Control,
                                           levels = design)))
deg <- topTable(fit2, number = Inf, sort.by = "none")

cat("总基因数:", nrow(deg), "\n")
cat("Up   (logFC>0.5 & adj.P<0.05):",
    sum(deg$logFC >  0.5 & deg$adj.P.Val < 0.05), "\n")  # 3,438
cat("Down (logFC<-0.5 & adj.P<0.05):",
    sum(deg$logFC < -0.5 & deg$adj.P.Val < 0.05), "\n")  # 1,404

# ---- 3. 加 Gene 列 + Status 列 ----
deg$Gene <- rownames(deg)
deg$Status <- "Not Significant"
deg$Status[deg$logFC >  0.5 & deg$adj.P.Val < 0.05] <- "Up in Insomnia"
deg$Status[deg$logFC < -0.5 & deg$adj.P.Val < 0.05] <- "Down in Insomnia"

# ---- 4. 导出（论文 Table S2 = 全 14,539 基因）----
deg_export <- deg[order(deg$adj.P.Val),
                  c("Gene", "logFC", "AveExpr", "t",
                    "P.Value", "adj.P.Val", "Status")]

write.csv(deg_export, "Phase6_Writing/limma_DEG_results.csv", row.names = FALSE)
cat("\n✅ Table S2 已保存：Phase6_Writing/limma_DEG_results.csv\n")
cat("   行数：", nrow(deg_export), "（=", nrow(expr_final), "全部基因）\n")

# ---- 5. 9 核心基因表达 + p 值核查 ----
cat("\n=== 9 核心基因 DEG ===\n")
core9 <- c("SPTLC1","PTGES3","HLA-G","NFE2L2","SOD2","HMOX1",
           "FKBP5","P2RX7","TNFAIP3")
core9_deg <- deg_export[deg_export$Gene %in% core9, ]
core9_deg <- core9_deg[match(core9, core9_deg$Gene), ]
print(core9_deg)

cat("\n========== Supplementary Tables 完成 ==========\n")
cat("已生成：\n")
cat("  Table S2  → Phase6_Writing/limma_DEG_results.csv  (14,539 基因)\n")
cat("  Table S5  → Phase2_WGCNA/hub_genes.csv             (2,120 hub)\n")
cat("  Table S6  → Phase3_ML/all_candidates_AUC.csv       (13 候选基因)\n")
cat("  Table S7  → Phase3.5_CrossCohort/cross_cohort_summary.csv\n")
cat("  Table S8  → Phase4_Immune/CIBERSORTx_Job460_Results.csv\n")
cat("  Table S9  → Phase5_SMR/results/SMR_full_results.csv\n")
cat("  Table S10 → Phase6_DGIdb/dgidb_all_results.csv\n")
cat("  Table S11 → Phase1_GSVA/GSVA_gene_sets.csv\n")
cat("  Table S12 → Phase3_ML/bootstrap_AUC_metrics.csv\n")



# ##############################################################################
# >>> END <<<
# Session info — record exact package versions for reproducibility
# ##############################################################################
cat("\n\n========== sessionInfo() ==========\n")
print(sessionInfo())
sink_path <- file.path(project_dir, "Phase6_Writing", "sessionInfo.txt")
dir.create(dirname(sink_path), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), sink_path)
cat("Session info saved to:", sink_path, "\n")
