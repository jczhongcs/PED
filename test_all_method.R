# ==============================================================================
# SCoPE2 核心补充实验：5 种策略对比 (修复样本对齐版)
# ==============================================================================

library(dplyr)
library(matrixStats)
library(impute)         # k-NN 插补
library(limma)          # 线性模型
library(preprocessCore) # 分位数归一化
library(batchelor)      # fastMNN
library(SingleCellExperiment)

# ------------------------------------------------------------------------------
# 1. 数据读取 (使用弹框选择)
# ------------------------------------------------------------------------------

# --- 选择 Cells.csv ---
message(">>> 请在弹出的窗口中选择: Cells.csv")
cells_path <- file.choose()
message(paste("已选择:", cells_path))
cells <- read.csv(cells_path, row.names = 1)

# --- 选择 Peptides-raw.csv ---
message(">>> 请在弹出的窗口中选择: Peptides-raw.csv")
peps_path <- file.choose()
message(paste("已选择:", peps_path))
peptides <- read.csv(peps_path)

# ------------------------------------------------------------------------------
# 2. 智能预处理 (自动修复方向)
# ------------------------------------------------------------------------------
message("\n>>> 正在预处理...")

# [关键修复] 检查 Cells.csv 是否需要转置
# 这里的逻辑是：看肽段表里的样本名是在 Cells 的列里，还是行里
pep_cols <- colnames(peptides)
if (length(intersect(colnames(cells), pep_cols)) > length(intersect(rownames(cells), pep_cols))) {
  message("   - 检测到 Cells.csv 样本在列上，正在执行转置...")
  cells <- as.data.frame(t(cells))
}

# 执行聚合 (Peptide -> Protein)
message("   - 执行聚合 (中位数)...")
raw_mat <- peptides %>%
  select(-peptide) %>%
  group_by(protein) %>%
  summarise(across(everything(), median, na.rm = TRUE)) %>%
  as.data.frame()
rownames(raw_mat) <- raw_mat$protein
raw_mat$protein <- NULL

# 对齐样本 (Intersection)
common <- intersect(colnames(raw_mat), rownames(cells))
message(sprintf("   - 成功对齐样本数: %d", length(common)))

if (length(common) == 0) {
  # 打印一些调试信息帮助定位
  message("DEBUG: Cells 行名示例: ", paste(head(rownames(cells)), collapse=", "))
  message("DEBUG: Protein 列名示例: ", paste(head(colnames(raw_mat)), collapse=", "))
  stop("错误: 样本仍然无法对齐！请检查上述示例是否匹配。")
}

raw_mat <- as.matrix(raw_mat[, common])
meta <- cells[common, ]

# 过滤 (保留 >50% 有效值)
keep <- rowMeans(!is.na(raw_mat)) > 0.5
data_filt <- raw_mat[keep, ]
message(sprintf("   - 过滤后保留蛋白数: %d", nrow(data_filt)))

message(">>> 执行 k-NN 插补 (k=3)...")
data_imputed <- impute.knn(data_filt, k = 3)$data

# Log2 转换检查
if(min(data_imputed, na.rm=T) >= 0) {
  message("   - 检测到线性数据，执行 Log2 转换...")
  data_imputed <- log2(data_imputed + 1e-6)
}

# ------------------------------------------------------------------------------
# 3. 生成 5 种校正数据
# ------------------------------------------------------------------------------

# [1] Ratio Only (Baseline)
write.csv(data_imputed, "Bench_01_Ratio.csv")
message("✅ [1/5] Ratio (Baseline) 生成完毕")

# [2] Median Centering
message(">>> 运行 Median Centering...")
data_median <- data_imputed
batches <- unique(meta$batch_chromatography)
g_med <- median(data_imputed, na.rm=TRUE)
for(b in batches) {
  idx <- which(meta$batch_chromatography == b)
  b_med <- median(data_imputed[, idx], na.rm=TRUE)
  data_median[, idx] <- data_imputed[, idx] - b_med + g_med
}
write.csv(data_median, "Bench_02_Median.csv")
message("✅ [2/5] Median 生成完毕")

# [3] Quantile Normalization
message(">>> 运行 Quantile Normalization...")
data_quant <- normalize.quantiles(data_imputed)
rownames(data_quant) <- rownames(data_imputed)
colnames(data_quant) <- colnames(data_imputed)
write.csv(data_quant, "Bench_03_Quantile.csv")
message("✅ [3/5] Quantile 生成完毕")

# [4] Limma (Blind)
message(">>> 运行 Limma...")
data_limma <- removeBatchEffect(data_imputed, batch = meta$batch_chromatography)
write.csv(data_limma, "Bench_04_Limma.csv")
message("✅ [4/5] Limma 生成完毕")

# [5] fastMNN (单细胞互近邻)
message(">>> 运行 fastMNN...")
tryCatch({
  fmnn_out <- fastMNN(data_imputed, batch = as.factor(meta$batch_chromatography))
  data_fastmnn <- assay(fmnn_out, "reconstructed")
  write.csv(as.matrix(data_fastmnn), "Bench_05_fastMNN.csv")
  message("✅ [5/5] fastMNN 生成完毕")
}, error = function(e) {
  message("⚠️ fastMNN 运行出错 (不影响其他结果): ", e$message)
})

message("\n🎉 全部完成！请转到 Python 进行 PED 评分。")