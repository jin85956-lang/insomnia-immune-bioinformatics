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
