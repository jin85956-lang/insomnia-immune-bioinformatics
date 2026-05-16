# ============================================================
# RUN_ALL.R — 一键运行全部分析（按顺序）
# ============================================================
# 用法：在 RStudio 中 source 整个文件，或逐行运行
# 注意：
#   1. 部分步骤需要联网（GEO 下载、eQTLGen、KEGG）
#   2. Phase 4 (CIBERSORTx) 和 Phase 5 (SMR) 需要中间手动操作
#   3. 完整流程在普通笔记本上约需 4-8 小时
# ============================================================

# 设置项目根目录
project_dir <- "F:/XiaoChen_Course/Insomnia_Project"
code_dir    <- file.path(project_dir, "code")   # 假设代码放在 code/ 下

# 按顺序运行
scripts <- c(
  "00_setup_and_data.R",        # 数据下载 + 预处理 (10 min)
  "01_Phase1_GSVA.R",           # GSVA 通路 (5 min)
  "02_Phase2_WGCNA.R",          # WGCNA (15-30 min)
  "03_Phase3_ML.R",             # ML LASSO/RF/Boruta (10 min)
  "04_Phase3.5_CrossCohort.R",  # 跨队列 (30 min, 含下载)
  "05_Phase0_eQTL.R",           # eQTLGen (15 min, 含下载)
  "06_Phase4_Immune.R",         # 阶段 A: CIBERSORTx 输入（之后需手动跑网站）
                                # 阶段 B: 结果解析（CIBERSORTx 跑完后再运行）
  "07_Phase5_SMR.R",            # SMR (需先用 SMR CLI 工具跑完)
  "08_Phase6_DGIdb.R",          # DGIdb（需先到网站查询）
  "09_Enrichment_GO_KEGG.R",    # GO/KEGG (10 min)
  "10_Bootstrap_AUC.R",         # Bootstrap (10 min)
  "11_Supplementary_Tables.R"   # 补充表
)

cat("========================================\n")
cat("  完整分析流程（共", length(scripts), "步）\n")
cat("========================================\n\n")

for (i in seq_along(scripts)) {
  cat(sprintf("[%d/%d] %s\n", i, length(scripts), scripts[i]))
  cat("----------------------------------------\n")
  source(file.path(code_dir, scripts[i]))
  cat("\n\n")
}

cat("========================================\n")
cat("  ✅ 全部分析完成\n")
cat("========================================\n")
cat("下一步：运行图片代码（你已有的 Untitled10）\n")
