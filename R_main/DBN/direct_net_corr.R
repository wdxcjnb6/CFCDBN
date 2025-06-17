library(R.matlab)
library(ggplot2)

# ==== 1. 读取评分和排序索引 ====
score_mat <- readMat("E:/PAC_network/Result/Direct_net_result/Direct_net_figure/Data/score_v6.mat")
v <- as.vector(score_mat$score.vals)

sorted_mat <- readMat("E:/PAC_network/Result/Direct_net_result/Direct_net_figure/Data/R_data_v6.mat")
sorted_indices <- as.vector(sorted_mat$sorted.indices)

# ==== 2. 设置路径 ====
source_net_path <- "E:/PAC_network/Result/Direct_net_result/Direct_net_data/TE/Anxiety"
mat_files <- list.files(source_net_path, pattern = "*_te.mat", full.names = TRUE)
num_subjects <- length(mat_files)

# ==== 3. 读取第一个样本，推断脑区数量 ====
example <- readMat(mat_files[1])
data_vec <- example$te.mat[[2,1,1]]
n_regions <- sqrt(length(data_vec))
if (n_regions %% 1 != 0) stop("数据长度无法开平方，可能不是方阵")
n_regions <- as.integer(n_regions)

# ==== 4. 构建网络张量 ====
all_networks <- array(NA, dim = c(n_regions, n_regions, num_subjects))
valid_count <- 0

for (i in seq_along(mat_files)) {
  subject <- readMat(mat_files[i])
  if (is.null(subject$te.mat)) {
    cat(sprintf("⚠️  %s 结构异常，跳过\n", basename(mat_files[i])))
    next
  }
  
  data_vec <- subject$te.mat[[2,1,1]]
  if (is.null(data_vec)) {
    cat(sprintf("⚠️  %s 中 data 为空，跳过\n", basename(mat_files[i])))
    next
  }
  
  data_mat <- matrix(data_vec, nrow = n_regions, ncol = n_regions)
  data_sorted <- data_mat[sorted_indices, sorted_indices]
  
  valid_count <- valid_count + 1
  all_networks[,,valid_count] <- data_sorted
}

# 裁剪有效被试
if (valid_count < num_subjects) {
  all_networks <- all_networks[, , 1:valid_count]
  cat(sprintf("✅ 有效被试数：%d/%d\n", valid_count, num_subjects))
}

# ==== 5. 逐条边与评分 Spearman 相关分析 ====
r_map <- matrix(NA, n_regions, n_regions)
p_map <- matrix(NA, n_regions, n_regions)

for (i in 1:n_regions) {
  for (j in 1:n_regions) {
    edge_vals <- all_networks[i, j, ]
    valid <- !is.na(edge_vals) & !is.na(v)
    if (sum(valid) >= 3) {
      test <- cor.test(edge_vals[valid], v[valid], method = "spearman")
      r_map[i, j] <- test$estimate
      p_map[i, j] <- test$p.value
    }
  }
}

# ==== 6. 提取显著连接 ====
sig_indices <- which(p_map < 0.05, arr.ind = TRUE)
num_sig <- nrow(sig_indices)
cat("✅ 共找到", num_sig, "个显著连接\n")

# ==== 7. 可视化每个显著连接 ====
output_dir <- "E:/PAC_network/R_code/有向脑网络/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (k in 1:num_sig) {
  i <- sig_indices[k, 1]
  j <- sig_indices[k, 2]
  edge_vals <- all_networks[i, j, ]
  
  df <- data.frame(score = v, edge = edge_vals)
  df$zero <- df$edge == 0
  non_zero <- df[!df$zero, ]
  
  p <- ggplot(df, aes(x = score, y = edge)) +
    geom_point(aes(color = zero), size = 2) +
    scale_color_manual(values = c("FALSE" = "#3399CC", "TRUE" = "#BBBBBB")) +
    geom_smooth(data = non_zero, method = "lm", se = FALSE, color = "red", linewidth = 1) +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("Edge (%d → %d) vs Score", i, j),
      subtitle = sprintf("r = %.3f, p = %.4g", r_map[i, j], p_map[i, j]),
      x = "Score", y = "Edge Strength"
    ) +
    theme(legend.position = "none")
  
  filename <- sprintf("Edge_%02d_to_%02d_p%.4f", i, j, p_map[i, j])
  ggsave(file.path(output_dir, paste0(filename, ".png")), p, width = 7, height = 6, dpi = 300)
}

# ==== 8. 保存结果 ====
save(r_map, p_map, file = file.path(output_dir, "r_p_map.RData"))
write.csv(r_map, file.path(output_dir, "r_map.csv"), row.names = FALSE)
write.csv(p_map, file.path(output_dir, "p_map.csv"), row.names = FALSE)

# ==== 显著边汇总表 ====
sig_table <- data.frame(
  Source = sig_indices[,1],
  Target = sig_indices[,2],
  Spearman_r = r_map[sig_indices],
  p_value = p_map[sig_indices]
)
write.csv(sig_table, file.path(output_dir, "SignificantEdges.csv"), row.names = FALSE)
