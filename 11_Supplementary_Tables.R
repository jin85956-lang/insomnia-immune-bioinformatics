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
