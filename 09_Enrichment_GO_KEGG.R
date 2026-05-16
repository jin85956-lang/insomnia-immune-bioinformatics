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
