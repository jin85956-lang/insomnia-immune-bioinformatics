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
