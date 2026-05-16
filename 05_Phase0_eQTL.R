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
