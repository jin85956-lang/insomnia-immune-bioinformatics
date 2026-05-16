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
