# ==== 加载必要包 ====
library(R.matlab)
library(pheatmap)
library(RColorBrewer)

# ==== 0. 配置路径 ====
mat_path <- "E:/PAC_network/Result/Direct_net_result/Direct_net_figure/Data/R_data.mat"
output_dir <- "E:/PAC_network/R_code/规模网络/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==== 1. 读取 .mat 数据 ====
mat <- readMat(mat_path)
sorted_indices <- as.vector(unlist(mat$sorted.indices))
group_net_data <- list(
  mat$group.net.data[[1]][[1]],
  mat$group.net.data[[2]][[1]]
)

# ==== 2. 功能网络信息 ====
network_names <- c("DMN", "SN", "FPN", "DAN", "SMN")
network_sizes <- c(10, 8, 6, 4, 6)
num_networks <- length(network_names)

roi_network_map <- c(rep("DMN",10),
                     rep("SN",8),
                     rep("FPN",6),
                     rep("DAN",4),
                     rep("SMN",6))

roi_network_sorted <- roi_network_map[sorted_indices]
network_indices <- lapply(network_names, function(net) which(roi_network_sorted == net))

# ==== 3. 初始化结果矩阵 ====
group1_mean_matrix <- matrix(NA, num_networks, num_networks)
group2_mean_matrix <- matrix(NA, num_networks, num_networks)
p_matrix <- matrix(NA, num_networks, num_networks)
t_matrix <- matrix(NA, num_networks, num_networks)

# ==== 4. 统计检验 ====
for (src in 1:num_networks) {
  for (tgt in 1:num_networks) {
    src_idx <- network_indices[[src]]
    tgt_idx <- network_indices[[tgt]]
    
    if (length(src_idx) == 0 || length(tgt_idx) == 0) next
    
    group1_vals <- apply(group_net_data[[1]][src_idx, tgt_idx, , drop = FALSE], 3, mean, na.rm = TRUE)
    group2_vals <- apply(group_net_data[[2]][src_idx, tgt_idx, , drop = FALSE], 3, mean, na.rm = TRUE)
    
    group1_mean_matrix[src, tgt] <- mean(group1_vals)
    group2_mean_matrix[src, tgt] <- mean(group2_vals)
    
    if (sum(!is.na(group1_vals)) >= 3 && sum(!is.na(group2_vals)) >= 3) {
      ttest <- t.test(group1_vals, group2_vals)
      p_matrix[src, tgt] <- ttest$p.value
      t_matrix[src, tgt] <- ttest$statistic
    }
  }
}

colnames(group1_mean_matrix) <- rownames(group1_mean_matrix) <- network_names
colnames(group2_mean_matrix) <- rownames(group2_mean_matrix) <- network_names
colnames(p_matrix) <- rownames(p_matrix) <- network_names

# ==== 5. 可视化参数设置 ====
purples <- colorRampPalette(brewer.pal(9, "Purples"))(100)
binary_colors <- c("white", "purple")
binary_breaks <- c(-0.1, 0.5, 1.1)

# ==== 6. 热图保存函数 ====
save_heatmap <- function(mat, file_base, title, color, breaks = NA) {
  # PNG
  png(file.path(output_dir, paste0(file_base, ".png")), width = 1200, height = 1000, res = 150)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
           main = title, color = color, breaks = breaks)
  dev.off()
  
  # SVG
  svg(file.path(output_dir, paste0(file_base, ".svg")), width = 8, height = 7)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
           main = title, color = color, breaks = breaks)
  dev.off()
}

# ==== 7. 生成并保存热图 ====

# Group1 & Group2
save_heatmap(group1_mean_matrix, "Group1_MeanConnectivity", "Group 1 Mean Connectivity", purples)
save_heatmap(group2_mean_matrix, "Group2_MeanConnectivity", "Group 2 Mean Connectivity", purples)

# p < 0.05 显著性图
sig_mat_raw <- (p_matrix < 0.05) * 1
save_heatmap(sig_mat_raw, "Significance_p_less_0.05", "Significant Connections (p < 0.05)", binary_colors, binary_breaks)

# ==== 8. FDR 校正 ====
p_vector <- as.vector(p_matrix)
valid_idx <- which(!is.na(p_vector))
p_adj <- rep(NA, length(p_vector))
if (length(valid_idx) > 0) {
  p_adj[valid_idx] <- p.adjust(p_vector[valid_idx], method = "fdr")
}
p_matrix_adj <- matrix(p_adj, nrow = num_networks, ncol = num_networks)
sig_matrix_adj_numeric <- (p_matrix_adj < 0.05) * 1

save_heatmap(sig_matrix_adj_numeric, "Significance_FDR_corrected", "FDR-corrected Significant Connections", binary_colors, binary_breaks)

# ==== 9. 输出统计结果表 ====
network_pairs <- as.vector(outer(network_names, network_names, function(a,b) paste(a, "→", b)))
network_stats <- data.frame(
  NetworkPair = network_pairs,
  Group1_Mean = as.vector(group1_mean_matrix),
  Group2_Mean = as.vector(group2_mean_matrix),
  p_value = as.vector(p_matrix),
  t_stat = as.vector(t_matrix),
  FDR_p = as.vector(p_matrix_adj)
)

write.csv(network_stats,
          file.path(output_dir, "组间显著性网络统计表.csv"),
          row.names = FALSE)

cat("原始 p < 0.05 的显著连接：\n")
print(subset(network_stats, !is.na(p_value) & p_value < 0.05))

cat("\nFDR 校正后 p < 0.05 的显著连接：\n")
print(subset(network_stats, !is.na(FDR_p) & FDR_p < 0.05))
