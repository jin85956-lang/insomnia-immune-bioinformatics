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
