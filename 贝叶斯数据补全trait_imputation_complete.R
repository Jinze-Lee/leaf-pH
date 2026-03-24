################################################################################
#  植物性状数据补全完整 R 实现
#  基于:
#    Joswig et al. (2023) Global Ecology and Biogeography 32:1395-1408
#    Schrodt et al. (2015) Global Ecology and Biogeography 24:1510-1521
#  方法: BHPMF + mice + missForest + 多指标评估
################################################################################


# ==============================================================================
# 第一部分: 环境准备与包安装
# ==============================================================================

# 安装所需包 (首次运行时取消注释)
# install.packages(c("mice", "missForest", "VIM", "naniar", "ggplot2",
#                    "reshape2", "vegan", "cluster", "corrplot",
#                    "dplyr", "tidyr", "purrr", "gridExtra"))
# install.packages("BHPMF")   # BHPMF 可能需要从 CRAN 或 GitHub 安装

# 加载包
library(BHPMF)       # 核心补全包: Bayesian Hierarchical Probabilistic Matrix Factorization
library(mice)        # 多重插补: Multivariate Imputation by Chained Equations
library(missForest)  # 随机森林插补
library(VIM)         # 缺失值可视化
library(naniar)      # 缺失值分析
library(ggplot2)     # 绘图
library(reshape2)    # 数据整形
library(vegan)       # Procrustes 检验
library(cluster)     # 轮廓系数分析
library(corrplot)    # 相关矩阵可视化
library(dplyr)       # 数据操作
library(tidyr)       # 数据整理
library(gridExtra)   # 多图排列


# ==============================================================================
# 第二部分: 数据载入与结构检查
# ==============================================================================

# ---- 2.1 使用 BHPMF 包内置示例数据 ----
data(trait.info)       # 植物性状矩阵 (行=个体, 列=性状), 缺失值为 NA
data(hierarchy.info)   # 分类层级矩阵 (列=层级: 门>纲>目>科>属>种>个体)

# 使用自有数据时, 替换为:
# trait.info      <- read.csv("your_trait_data.csv",     row.names = 1)
# hierarchy.info  <- read.csv("your_hierarchy_data.csv", row.names = 1)
# trait.info      <- as.matrix(trait.info)
# hierarchy.info  <- as.matrix(hierarchy.info)

cat("=== 数据基本信息 ===\n")
cat("性状矩阵维度:", dim(trait.info), "\n")
cat("层级矩阵维度:", dim(hierarchy.info), "\n")
cat("性状名称:", colnames(trait.info), "\n")
cat("层级列名:", colnames(hierarchy.info), "\n")

# ---- 2.2 数据类型确认 ----
trait.matrix     <- as.matrix(trait.info)
hierarchy.matrix <- as.matrix(hierarchy.info)

# 确保所有性状列均为数值型
storage.mode(trait.matrix) <- "double"


# ==============================================================================
# 第三部分: 缺失数据诊断
# ==============================================================================

# ---- 3.1 基本缺失统计 ----
cat("\n=== 缺失数据统计 ===\n")

miss_summary <- data.frame(
  trait        = colnames(trait.matrix),
  n_missing    = colSums(is.na(trait.matrix)),
  pct_missing  = round(colMeans(is.na(trait.matrix)) * 100, 2),
  n_observed   = colSums(!is.na(trait.matrix))
)
print(miss_summary)

overall_miss <- mean(is.na(trait.matrix)) * 100
cat(sprintf("\n整体缺失率: %.2f%%\n", overall_miss))

# ---- 3.2 每行缺失情况 ----
row_miss <- rowSums(is.na(trait.matrix))
cat(sprintf("有完整观测的行: %d (%.1f%%)\n",
            sum(row_miss == 0),
            sum(row_miss == 0) / nrow(trait.matrix) * 100))
cat(sprintf("全部缺失的行:   %d\n", sum(row_miss == ncol(trait.matrix))))

# ---- 3.3 可视化缺失模式 ----
# 方法1: 使用 naniar 包
p_miss_upset <- gg_miss_upset(as.data.frame(trait.matrix), nsets = 10)
print(p_miss_upset)

# 方法2: 使用 VIM 包的聚合图
aggr_plot <- aggr(as.data.frame(trait.matrix),
                  col   = c("navyblue", "red"),
                  numbers = TRUE,
                  sortVars = TRUE,
                  labels = colnames(trait.matrix),
                  cex.axis = 0.7,
                  gap = 3,
                  ylab = c("缺失值直方图", "缺失模式"))

# 方法3: 热图展示缺失模式
miss_pattern <- is.na(trait.matrix)
miss_df      <- melt(miss_pattern)
colnames(miss_df) <- c("Observation", "Trait", "Missing")

p_heatmap <- ggplot(miss_df, aes(x = Trait, y = Observation, fill = Missing)) +
  geom_tile() +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                    labels = c("已观测", "缺失")) +
  labs(title = "缺失值热图", x = "性状", y = "观测个体", fill = "状态") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
print(p_heatmap)

# ---- 3.4 缺失机制检验 (MCAR/MAR/MNAR) ----
# Little's MCAR 检验 (需要 naniar 包中的 mcar_test)
# 原假设: 数据完全随机缺失 (MCAR)
if (requireNamespace("naniar", quietly = TRUE)) {
  mcar_result <- tryCatch(
    mcar_test(as.data.frame(trait.matrix)),
    error = function(e) NULL
  )
  if (!is.null(mcar_result)) {
    cat("\n=== Little's MCAR 检验 ===\n")
    print(mcar_result)
    cat("若 p > 0.05: 数据可能为 MCAR (完全随机缺失)\n")
    cat("若 p < 0.05: 数据可能为 MAR 或 MNAR\n")
  }
}


# ==============================================================================
# 第四部分: 数据变换
# (依据 Joswig et al. 2023 及 BHPMF 包文档: 先取 log, 再 z 变换)
# ==============================================================================

# ---- 4.1 确认所有性状为正值 (log 变换前提) ----
neg_check <- apply(trait.matrix, 2, function(x) any(x <= 0, na.rm = TRUE))
if (any(neg_check)) {
  warning("以下列含非正值, log 变换前需处理: ",
          paste(colnames(trait.matrix)[neg_check], collapse = ", "))
}

# ---- 4.2 log 变换 ----
trait_log <- log(trait.matrix)   # 自然对数; 若数据含0, 改用 log(x + 1)

# ---- 4.3 z 变换 (按列标准化) ----
# z_log(y) = (log(y) - mu(log(k))) / sigma(log(k))
trait_mean_log <- colMeans(trait_log, na.rm = TRUE)
trait_sd_log   <- apply(trait_log, 2, sd, na.rm = TRUE)

trait_zlog <- sweep(trait_log, 2, trait_mean_log, "-")
trait_zlog <- sweep(trait_zlog, 2, trait_sd_log,  "/")

cat("\n=== 变换后各性状统计 (z-log 空间) ===\n")
print(round(apply(trait_zlog, 2, function(x)
  c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE),
    min  = min(x,  na.rm=TRUE), max = max(x, na.rm=TRUE))), 3))

# ---- 4.4 反变换函数 (补全后恢复原尺度) ----
backtransform <- function(zlog_mat, mean_log, sd_log) {
  # 反 z 变换
  log_mat <- sweep(zlog_mat, 2, sd_log,   "*")
  log_mat <- sweep(log_mat,  2, mean_log, "+")
  # 反 log 变换
  exp(log_mat)
}


# ==============================================================================
# 第五部分: 数据打孔 (Perforation) — 用于评估补全精度
# (依据 Joswig et al. 2023: 随机删除已有观测值, 保留至少每行/列各一个)
# ==============================================================================

perforate_data <- function(X, miss_rate = 0.5, seed = 42) {
  # 随机删除 miss_rate 比例的已观测值
  # 约束: 每行至少保留 1 个观测, 每列至少保留 1 个观测
  set.seed(seed)
  X_sparse <- X
  obs_idx  <- which(!is.na(X), arr.ind = TRUE)

  # 计算目标删除数量
  n_remove <- floor(nrow(obs_idx) * miss_rate)

  # 逐步随机删除, 确保行列约束
  removed <- 0
  shuffled <- sample(nrow(obs_idx))
  for (i in shuffled) {
    r <- obs_idx[i, 1]; c <- obs_idx[i, 2]
    row_obs <- sum(!is.na(X_sparse[r, ]))
    col_obs <- sum(!is.na(X_sparse[, c]))
    if (row_obs > 1 && col_obs > 1) {
      X_sparse[r, c] <- NA
      removed <- removed + 1
      if (removed >= n_remove) break
    }
  }
  actual_miss <- mean(is.na(X_sparse)) * 100
  cat(sprintf("打孔后缺失率: %.2f%%\n", actual_miss))
  return(X_sparse)
}

# 生成稀疏版本 (50% 额外缺失), 重复 3 次 (依据论文)
trait_sparse_list <- lapply(1:3, function(rep) {
  perforate_data(trait_zlog, miss_rate = 0.5, seed = rep * 100)
})
trait_sparse <- trait_sparse_list[[1]]  # 使用第一次重复作为主数据集


# ==============================================================================
# 第六部分: BHPMF 补全 (核心方法)
# (依据 Schrodt et al. 2015; Joswig et al. 2023)
# ==============================================================================

# ---- 6.1 输出路径设置 ----
output_dir <- tempdir()
mean_path  <- file.path(output_dir, "mean_gap_filled.txt")
std_path   <- file.path(output_dir, "std_gap_filled.txt")
tmp_dir    <- file.path(output_dir, "BHPMF_tmp")
dir.create(tmp_dir, showWarnings = FALSE)

# ---- 6.2 运行 BHPMF (Joswig et al. 参数: 1000次迭代, 前200丢弃, 每20次取一) ----
cat("\n=== 运行 BHPMF 补全 ===\n")
GapFilling(
  X                        = trait_sparse,          # 含缺失值的性状矩阵 (z-log 空间)
  hierarchy.info           = hierarchy.matrix,       # 分类层级矩阵
  prediction.level         = ncol(hierarchy.matrix), # 在最底层 (个体) 补全
  used.num.hierarchy.levels = ncol(hierarchy.matrix) - 1, # 使用所有层级信息
  num.samples              = 1000,  # Gibbs 采样总次数
  burn                     = 200,   # 前200次丢弃 (burn-in)
  gaps                     = 20,    # 每20次保留1次 (避免自相关) → 40个有效样本
  num.latent.feats         = 10,    # 潜在因子数量
  tuning                   = FALSE, # 不自动调参 (生产环境可设为 TRUE)
  tmp.dir                  = tmp_dir,
  mean.gap.filled.output.path = mean_path,
  std.gap.filled.output.path  = std_path,
  rmse.plot.test.data      = TRUE,
  verbose                  = TRUE
)

# ---- 6.3 读取 BHPMF 结果 ----
imputed_mean_zlog <- as.matrix(read.table(mean_path, header = FALSE,
                                          sep = "\t", na.strings = "NA"))
imputed_std_zlog  <- as.matrix(read.table(std_path,  header = FALSE,
                                          sep = "\t", na.strings = "NA"))

colnames(imputed_mean_zlog) <- colnames(trait.matrix)
colnames(imputed_std_zlog)  <- colnames(trait.matrix)

# ---- 6.4 反变换回原始尺度 ----
imputed_mean_orig <- backtransform(imputed_mean_zlog, trait_mean_log, trait_sd_log)

cat("\n=== BHPMF 补全完成 ===\n")
cat("补全矩阵维度:", dim(imputed_mean_zlog), "\n")
cat("残余 NA 数量:", sum(is.na(imputed_mean_zlog)), "\n")

# ---- 6.5 参数调优 (可选, 耗时) ----
# tune_result <- TuneBhpmf(
#   X                        = trait_sparse,
#   hierarchy.info           = hierarchy.matrix,
#   num.folds                = 5,
#   num.samples              = 500,
#   burn                     = 100,
#   gaps                     = 10,
#   tmp.dir                  = tmp_dir,
#   verbose                  = TRUE
# )
# best_latent_feats <- tune_result$best.num.latent.feats
# cat("最优潜在因子数:", best_latent_feats, "\n")

# ---- 6.6 交叉验证 RMSE (可选) ----
# cv_result <- CalculateCvRmse(
#   X                        = trait_sparse,
#   hierarchy.info           = hierarchy.matrix,
#   num.folds                = 10,
#   num.samples              = 500,
#   burn                     = 100,
#   gaps                     = 10,
#   num.latent.feats         = 10,
#   verbose                  = TRUE
# )
# cat(sprintf("CV RMSE: %.4f ± %.4f\n", cv_result$avg.rmse, cv_result$std.rmse))


# ==============================================================================
# 第七部分: 其他补全方法 (对比基准)
# ==============================================================================

# ---- 7.1 mice (多重链式方程插补) ----
cat("\n=== 运行 mice 补全 ===\n")
mice_result <- mice(
  data     = as.data.frame(trait_sparse),
  m        = 5,          # 生成 5 个插补数据集
  method   = "pmm",      # 预测均值匹配 (适合连续变量)
  maxit    = 50,          # 最大迭代次数
  seed     = 123,
  printFlag = FALSE
)
# 提取第一个完整数据集
imputed_mice <- complete(mice_result, action = 1)
imputed_mice <- as.matrix(imputed_mice)

# 汇总多重插补结果 (取5个数据集的均值)
imputed_mice_mean <- Reduce("+", lapply(1:5, function(i)
  as.matrix(complete(mice_result, action = i)))) / 5

cat("mice 补全完成, 残余 NA:", sum(is.na(imputed_mice)), "\n")

# ---- 7.2 missForest (随机森林插补) ----
cat("\n=== 运行 missForest 补全 ===\n")
mf_result <- missForest(
  xmis     = trait_sparse,
  maxiter  = 10,
  ntree    = 100,
  verbose  = TRUE
)
imputed_mf <- mf_result$ximp
cat(sprintf("missForest OOB 误差: %.4f\n", mf_result$OOBerror[1]))

# ---- 7.3 均值插补 (最简单基准) ----
imputed_mean_simple <- trait_sparse
for (j in seq_len(ncol(trait_sparse))) {
  col_mean <- mean(trait_sparse[, j], na.rm = TRUE)
  imputed_mean_simple[is.na(trait_sparse[, j]), j] <- col_mean
}

# ---- 7.4 物种均值插补 (分类学基准) ----
# 需要 hierarchy.matrix 中包含物种列信息
# 此处以最后一列作为物种标识符
species_col <- hierarchy.matrix[, ncol(hierarchy.matrix)]
imputed_species_mean <- trait_sparse

for (j in seq_len(ncol(trait_sparse))) {
  miss_idx <- which(is.na(trait_sparse[, j]))
  for (i in miss_idx) {
    sp    <- species_col[i]
    sp_idx <- which(species_col == sp & !is.na(trait_sparse[, j]))
    if (length(sp_idx) > 0) {
      imputed_species_mean[i, j] <- mean(trait_sparse[sp_idx, j], na.rm = TRUE)
    } else {
      imputed_species_mean[i, j] <- mean(trait_sparse[, j], na.rm = TRUE)
    }
  }
}


# ==============================================================================
# 第八部分: 精度评估
# (依据 Joswig et al. 2023: RMSE、残差分布、性状相关、PCA 比较)
# ==============================================================================

# ---- 8.1 辅助函数 ----
calc_rmse <- function(observed, predicted, mask_na) {
  # mask_na: 指示哪些位置是被人为打孔的位置
  obs  <- observed[mask_na]
  pred <- predicted[mask_na]
  valid <- !is.na(obs) & !is.na(pred)
  sqrt(mean((obs[valid] - pred[valid])^2))
}

calc_rmse_per_trait <- function(observed, predicted, mask_na) {
  sapply(seq_len(ncol(observed)), function(j) {
    obs  <- observed[mask_na[, j], j]
    pred <- predicted[mask_na[, j], j]
    valid <- !is.na(obs) & !is.na(pred)
    if (sum(valid) == 0) return(NA)
    sqrt(mean((obs[valid] - pred[valid])^2))
  })
}

# ---- 8.2 计算打孔位置掩码 ----
# 原始有值但被打孔的位置
hole_mask <- !is.na(trait_zlog) & is.na(trait_sparse)
cat(sprintf("\n评估用打孔数量: %d\n", sum(hole_mask)))

# ---- 8.3 各方法 RMSE 比较 ----
methods_list <- list(
  BHPMF          = imputed_mean_zlog,
  mice           = imputed_mice_mean,
  missForest     = imputed_mf,
  均值插补       = imputed_mean_simple,
  物种均值插补   = imputed_species_mean
)

rmse_results <- data.frame(
  方法     = names(methods_list),
  总体RMSE = sapply(methods_list, function(imp)
    calc_rmse(trait_zlog, imp, hole_mask))
)
rmse_results <- rmse_results[order(rmse_results$总体RMSE), ]
cat("\n=== 各方法总体 RMSE (z-log 空间) ===\n")
print(rmse_results)

# 各性状 RMSE
rmse_per_trait <- sapply(methods_list, function(imp)
  calc_rmse_per_trait(trait_zlog, imp, hole_mask))
rownames(rmse_per_trait) <- colnames(trait_zlog)
cat("\n=== 各性状 RMSE ===\n")
print(round(rmse_per_trait, 4))

# ---- 8.4 残差分析 (仅 BHPMF) ----
# y_residual = y_imputed - y_observed (反变换后的原始尺度)
obs_orig  <- backtransform(trait_zlog,         trait_mean_log, trait_sd_log)
imp_orig  <- imputed_mean_orig

residuals_orig <- imp_orig - obs_orig  # 只关心打孔位置
residuals_vec  <- residuals_orig[hole_mask]

cat("\n=== BHPMF 残差统计 (原始尺度) ===\n")
cat(sprintf("均值:   %.4f\n", mean(residuals_vec, na.rm=TRUE)))
cat(sprintf("标准差: %.4f\n", sd(residuals_vec,   na.rm=TRUE)))
cat(sprintf("中位数: %.4f\n", median(residuals_vec, na.rm=TRUE)))

# 残差图
p_resid <- ggplot(data.frame(residuals = residuals_vec), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 50,
                 fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "BHPMF 补全残差分布 (原始尺度)",
       x = "残差 (补全值 - 观测值)", y = "密度") +
  theme_minimal()
print(p_resid)

# 不确定性 vs RMSE 散点图
uncert_vec   <- imputed_std_zlog[hole_mask]
abs_err_vec  <- abs((imputed_mean_zlog - trait_zlog)[hole_mask])
valid_idx    <- !is.na(uncert_vec) & !is.na(abs_err_vec)

p_uncert <- ggplot(data.frame(sd = uncert_vec[valid_idx],
                               abs_err = abs_err_vec[valid_idx]),
                   aes(x = sd, y = abs_err)) +
  geom_point(alpha = 0.3, size = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "BHPMF 不确定性 (SD) vs 绝对误差",
       x = "补全标准差 (SD)", y = "绝对误差") +
  theme_minimal()
print(p_uncert)


# ==============================================================================
# 第九部分: 性状-性状相关分析
# (依据 Joswig et al. 2023: Pearson 相关 + PCA + Procrustes 检验)
# ==============================================================================

# ---- 9.1 Pearson 相关矩阵比较 ----
cor_obs  <- cor(trait_zlog,         use = "pairwise.complete.obs", method = "pearson")
cor_imp  <- cor(imputed_mean_zlog,  use = "pairwise.complete.obs", method = "pearson")
cor_mice <- cor(imputed_mice_mean,  use = "pairwise.complete.obs", method = "pearson")

# 可视化相关矩阵
par(mfrow = c(1, 3))
corrplot(cor_obs,  method = "color", title = "观测数据",   mar = c(0,0,1,0))
corrplot(cor_imp,  method = "color", title = "BHPMF 补全", mar = c(0,0,1,0))
corrplot(cor_mice, method = "color", title = "mice 补全",   mar = c(0,0,1,0))
par(mfrow = c(1, 1))

# Pearson 相关系数均值比较
cat("\n=== Pearson 相关系数均值 ===\n")
upper_tri <- upper.tri(cor_obs)
cat(sprintf("观测数据:   %.4f\n", mean(abs(cor_obs[upper_tri]),  na.rm=TRUE)))
cat(sprintf("BHPMF 补全: %.4f\n", mean(abs(cor_imp[upper_tri]),  na.rm=TRUE)))
cat(sprintf("mice 补全:  %.4f\n", mean(abs(cor_mice[upper_tri]), na.rm=TRUE)))

# ---- 9.2 PCA 分析 ----
# 使用无缺失的完整补全矩阵做 PCA
pca_obs <- princomp(na.omit(trait_zlog))
pca_imp <- princomp(imputed_mean_zlog[complete.cases(imputed_mean_zlog), ])

cat("\n=== 观测数据 PCA 解释方差 ===\n")
print(summary(pca_obs))

# ---- 9.3 Procrustes 检验 (依据 Joswig et al. 2023) ----
# 使用共同行进行比较
common_rows <- which(complete.cases(trait_zlog) & complete.cases(imputed_mean_zlog))
if (length(common_rows) > 10) {
  pca_obs_sub <- prcomp(trait_zlog[common_rows, ], scale. = FALSE)
  pca_imp_sub <- prcomp(imputed_mean_zlog[common_rows, ], scale. = FALSE)

  # Procrustes 分析: 旋转使两者最相似
  proc_result  <- procrustes(pca_obs_sub$x[, 1:2], pca_imp_sub$x[, 1:2])
  proto_test   <- protest(pca_obs_sub$x[, 1:2],    pca_imp_sub$x[, 1:2],
                           permutations = 999)

  cat("\n=== Procrustes 检验结果 ===\n")
  cat(sprintf("m12 统计量 (残差): %.4f\n", proc_result$ss))
  cat(sprintf("相关性 (r):         %.4f\n", sqrt(1 - proc_result$ss)))
  cat(sprintf("显著性 p 值:        %.4f\n", proto_test$signif))
  cat("p < 0.05 表示两个 PCA 结果非随机相似\n")
}


# ==============================================================================
# 第十部分: 分类聚类分析
# (依据 Joswig et al. 2023: 轮廓系数、变异系数、到聚类均值的距离)
# ==============================================================================

# ---- 10.1 到物种均值的距离 ----
calc_dist_to_cluster_mean <- function(X, species_labels) {
  # 计算每个观测值到其物种均值的距离 (z-log 空间)
  species_means <- aggregate(X, by = list(species = species_labels), FUN = mean, na.rm = TRUE)
  rownames(species_means) <- species_means$species
  species_means$species   <- NULL

  distances <- matrix(NA, nrow(X), ncol(X))
  for (i in seq_len(nrow(X))) {
    sp <- species_labels[i]
    if (sp %in% rownames(species_means)) {
      distances[i, ] <- as.numeric(X[i, ]) - as.numeric(species_means[sp, ])
    }
  }
  return(distances)
}

# 计算观测数据和补全数据的到物种均值距离
species_labels <- as.character(species_col)
dist_obs <- calc_dist_to_cluster_mean(trait_zlog,        species_labels)
dist_imp <- calc_dist_to_cluster_mean(imputed_mean_zlog, species_labels)

# 比较分布
cat("\n=== 到物种均值的平均距离 (z-log 空间) ===\n")
cat(sprintf("观测数据:   %.4f\n", mean(abs(dist_obs), na.rm=TRUE)))
cat(sprintf("BHPMF 补全: %.4f\n", mean(abs(dist_imp), na.rm=TRUE)))

# ---- 10.2 变异系数 (CV) 分析 ----
# 依据 Joswig et al. 2023: CV = SD/mean, 基于反变换后的原始尺度
calc_cv_by_species <- function(X_orig, species_labels) {
  # X_orig: 原始尺度矩阵
  cv_list <- lapply(colnames(X_orig), function(trait) {
    df <- data.frame(val = X_orig[, trait], sp = species_labels)
    df <- df[!is.na(df$val), ]
    sp_cv <- tapply(df$val, df$sp, function(x) {
      if (length(x) > 1) sd(x) / mean(x) else NA
    })
    data.frame(trait = trait, cv = sp_cv[!is.na(sp_cv)])
  })
  do.call(rbind, cv_list)
}

cv_obs <- calc_cv_by_species(obs_orig,   species_labels)
cv_imp <- calc_cv_by_species(imp_orig,   species_labels)

cv_summary <- data.frame(
  trait    = unique(cv_obs$trait),
  cv_obs   = tapply(cv_obs$cv, cv_obs$trait, mean, na.rm=TRUE),
  cv_imp   = tapply(cv_imp$cv, cv_imp$trait, mean, na.rm=TRUE)
)
cv_summary$cv_ratio <- cv_summary$cv_imp / cv_summary$cv_obs

cat("\n=== 物种内变异系数 (CV) 比较 (原始尺度) ===\n")
print(round(cv_summary, 4))
cat("cv_ratio > 1 表示补全数据物种内变异增大 (聚类偏差)\n")

# ---- 10.3 轮廓系数 (多变量聚类质量) ----
if (length(common_rows) > 20) {
  # 仅使用完整观测行
  obs_complete <- trait_zlog[complete.cases(trait_zlog), ]
  imp_complete <- imputed_mean_zlog[complete.cases(imputed_mean_zlog), ]

  # 以物种标签作为聚类标签
  sp_complete_obs <- species_labels[complete.cases(trait_zlog)]
  sp_complete_imp <- species_labels[complete.cases(imputed_mean_zlog)]

  # 将物种标签转为整数
  sp_int_obs <- as.integer(factor(sp_complete_obs))
  sp_int_imp <- as.integer(factor(sp_complete_imp))

  # 计算轮廓系数 (使用子集避免计算过慢)
  set.seed(42)
  n_sub <- min(500, nrow(obs_complete))
  sub_idx <- sample(nrow(obs_complete), n_sub)

  dist_obs_sub <- dist(obs_complete[sub_idx, ])
  sil_obs <- silhouette(sp_int_obs[sub_idx], dist_obs_sub)

  dist_imp_sub <- dist(imp_complete[sub_idx, ])
  sil_imp <- silhouette(sp_int_imp[sub_idx], dist_imp_sub)

  cat(sprintf("\n物种轮廓系数 - 观测数据:   %.4f\n", mean(sil_obs[, 3])))
  cat(sprintf("物种轮廓系数 - BHPMF 补全: %.4f\n", mean(sil_imp[, 3])))
  cat("值越高表示物种内聚类越紧密; BHPMF 可能人为增大该值 (聚类偏差)\n")
}


# ==============================================================================
# 第十一部分: 性状分布可视化
# (依据 Joswig et al. 2023 图3: 观测 vs 补全的密度曲线)
# ==============================================================================

# ---- 11.1 密度比较图 ----
plot_trait_dist <- function(trait_idx, obs_orig, imp_orig, sparse_orig, trait_name) {
  obs_vals    <- obs_orig[!is.na(obs_orig[, trait_idx]), trait_idx]
  imp_vals    <- imp_orig[, trait_idx]
  sparse_vals <- sparse_orig[!is.na(sparse_orig[, trait_idx]), trait_idx]

  df_plot <- rbind(
    data.frame(val = obs_vals,    source = "观测 (OBS)"),
    data.frame(val = imp_vals,    source = "BHPMF 补全"),
    data.frame(val = sparse_vals, source = "稀疏观测")
  )

  ggplot(df_plot, aes(x = val, color = source, fill = source)) +
    geom_density(alpha = 0.2, linewidth = 0.8) +
    labs(title = trait_name, x = "性状值", y = "密度", color = "", fill = "") +
    scale_color_manual(values = c("观测 (OBS)" = "steelblue",
                                   "BHPMF 补全" = "darkorange",
                                   "稀疏观测"   = "darkblue")) +
    scale_fill_manual(values  = c("观测 (OBS)" = "steelblue",
                                   "BHPMF 补全" = "darkorange",
                                   "稀疏观测"   = "darkblue")) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# 将稀疏数据也反变换回原始尺度
sparse_orig <- backtransform(trait_sparse, trait_mean_log, trait_sd_log)

# 绘制所有性状的分布图
plot_list <- lapply(seq_len(ncol(trait.matrix)), function(j)
  plot_trait_dist(j, obs_orig, imp_orig, sparse_orig, colnames(trait.matrix)[j])
)
do.call(grid.arrange, c(plot_list, list(ncol = 3)))

# ---- 11.2 均值比较条形图 ----
means_df <- data.frame(
  trait  = colnames(trait.matrix),
  OBS    = colMeans(obs_orig,    na.rm=TRUE),
  BHPMF  = colMeans(imp_orig,    na.rm=TRUE),
  mice   = colMeans(backtransform(imputed_mice_mean, trait_mean_log, trait_sd_log), na.rm=TRUE)
)
means_long <- melt(means_df, id.vars = "trait", variable.name = "方法", value.name = "均值")

p_means <- ggplot(means_long, aes(x = trait, y = 均值, fill = 方法)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(title = "各方法性状均值比较 (原始尺度)",
       x = "性状", y = "均值") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_means)

# ---- 11.3 补全不确定性热图 ----
std_df <- as.data.frame(imputed_std_zlog)
std_df$row <- seq_len(nrow(std_df))
std_long <- melt(std_df, id.vars = "row", variable.name = "trait", value.name = "SD")

p_sd <- ggplot(std_long, aes(x = trait, y = row, fill = SD)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred",
                       name = "SD\n(不确定性)") +
  labs(title = "BHPMF 补全不确定性 (SD)", x = "性状", y = "观测个体") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
print(p_sd)


# ==============================================================================
# 第十二部分: 重复稳定性分析
# (依据 Joswig et al. 2023: 3 次重复评估稳定性)
# ==============================================================================

if (FALSE) {  # 耗时较长, 需要时设为 TRUE
  cat("\n=== 运行 3 次重复 BHPMF ===\n")
  rep_results <- list()
  for (rep in 1:3) {
    tmp_rep <- file.path(tempdir(), paste0("rep", rep))
    dir.create(tmp_rep, showWarnings = FALSE)
    mean_rep <- file.path(tmp_rep, "mean.txt")
    std_rep  <- file.path(tmp_rep, "std.txt")

    GapFilling(
      X                           = trait_sparse_list[[rep]],
      hierarchy.info              = hierarchy.matrix,
      num.samples                 = 1000,
      burn                        = 200,
      gaps                        = 20,
      num.latent.feats            = 10,
      tmp.dir                     = tmp_rep,
      mean.gap.filled.output.path = mean_rep,
      std.gap.filled.output.path  = std_rep,
      rmse.plot.test.data         = FALSE,
      verbose                     = FALSE
    )
    rep_results[[rep]] <- as.matrix(read.table(mean_rep, header=FALSE, sep="\t"))
  }

  # 计算各重复之间的相关性 (稳定性指标)
  cor_rep12 <- cor(as.vector(rep_results[[1]]), as.vector(rep_results[[2]]),
                   use="complete.obs")
  cor_rep13 <- cor(as.vector(rep_results[[1]]), as.vector(rep_results[[3]]),
                   use="complete.obs")
  cor_rep23 <- cor(as.vector(rep_results[[2]]), as.vector(rep_results[[3]]),
                   use="complete.obs")
  cat(sprintf("重复1-2 Pearson r: %.4f\n", cor_rep12))
  cat(sprintf("重复1-3 Pearson r: %.4f\n", cor_rep13))
  cat(sprintf("重复2-3 Pearson r: %.4f\n", cor_rep23))
}


# ==============================================================================
# 第十三部分: 不同缺失率下的 RMSE 曲线
# (依据 Joswig et al. 2023: 1-80% 缺失率评估)
# ==============================================================================

if (FALSE) {  # 耗时较长, 需要时设为 TRUE
  miss_rates   <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)
  rmse_by_rate <- numeric(length(miss_rates))

  for (k in seq_along(miss_rates)) {
    rate <- miss_rates[k]
    cat(sprintf("\n--- 缺失率 %.0f%% ---\n", rate * 100))

    tmp_k  <- file.path(tempdir(), paste0("miss_", k))
    dir.create(tmp_k, showWarnings = FALSE)
    X_k    <- perforate_data(trait_zlog, miss_rate = rate, seed = k)
    hole_k <- !is.na(trait_zlog) & is.na(X_k)
    m_path <- file.path(tmp_k, "mean.txt")
    s_path <- file.path(tmp_k, "std.txt")

    GapFilling(X = X_k, hierarchy.info = hierarchy.matrix,
               num.samples = 500, burn = 100, gaps = 10, num.latent.feats = 10,
               tmp.dir = tmp_k,
               mean.gap.filled.output.path = m_path,
               std.gap.filled.output.path  = s_path,
               rmse.plot.test.data = FALSE, verbose = FALSE)

    imp_k <- as.matrix(read.table(m_path, header=FALSE, sep="\t"))
    rmse_by_rate[k] <- calc_rmse(trait_zlog, imp_k, hole_k)
  }

  # 绘制 RMSE vs 缺失率曲线
  df_rmse_rate <- data.frame(miss_rate = miss_rates * 100, RMSE = rmse_by_rate)
  p_rmse_curve <- ggplot(df_rmse_rate, aes(x = miss_rate, y = RMSE)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "red", size = 3) +
    labs(title = "BHPMF RMSE 随缺失率变化",
         x = "缺失率 (%)", y = "RMSE (z-log 空间)") +
    theme_minimal()
  print(p_rmse_curve)
}


# ==============================================================================
# 第十四部分: 结果导出
# ==============================================================================

# ---- 14.1 保存 BHPMF 补全结果 (原始尺度) ----
write.csv(imputed_mean_orig,
          "BHPMF_imputed_mean_original_scale.csv", row.names = TRUE)
write.csv(imputed_mean_zlog,
          "BHPMF_imputed_mean_zlog_scale.csv",     row.names = TRUE)
write.csv(imputed_std_zlog,
          "BHPMF_imputed_std_zlog_scale.csv",      row.names = TRUE)

# ---- 14.2 保存 mice 补全结果 ----
write.csv(backtransform(imputed_mice_mean, trait_mean_log, trait_sd_log),
          "mice_imputed_original_scale.csv", row.names = TRUE)

# ---- 14.3 保存评估指标 ----
write.csv(rmse_results,    "rmse_comparison.csv",    row.names = FALSE)
write.csv(rmse_per_trait,  "rmse_per_trait.csv",     row.names = TRUE)
write.csv(cv_summary,      "cv_comparison.csv",      row.names = FALSE)

# ---- 14.4 综合报告 ----
cat("\n", strrep("=", 60), "\n")
cat("=== 补全完成! 结果摘要 ===\n")
cat(strrep("=", 60), "\n")
cat(sprintf("原始数据: %d 行 x %d 列\n",    nrow(trait.matrix), ncol(trait.matrix)))
cat(sprintf("整体缺失率:      %.2f%%\n",      overall_miss))
cat(sprintf("BHPMF 总体 RMSE: %.4f (z-log)\n", rmse_results$总体RMSE[rmse_results$方法=="BHPMF"]))
cat(sprintf("mice  总体 RMSE: %.4f (z-log)\n", rmse_results$总体RMSE[rmse_results$方法=="mice"]))
cat("\n推荐使用 BHPMF, 因为它:\n")
cat("  1. 利用分类层级结构信息\n")
cat("  2. 利用性状间相关结构\n")
cat("  3. 提供逐值的不确定性估计 (SD)\n")
cat("  4. 即使矩阵极稀疏也能稳健运行\n")
cat("\n注意偏差:\n")
cat("  - BHPMF 可能增大物种内聚类紧密度 (taxonomic clustering bias)\n")
cat("  - 建议对重要结论进行偏差检验 (见 Joswig et al. 2023)\n")
cat(strrep("=", 60), "\n")
