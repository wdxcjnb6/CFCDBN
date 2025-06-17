library(ggplot2)
library(readr)
library(stringr)

data_dir <- "E:/PAC_network/R_code/Corr"
plot_dir <- file.path(data_dir, "EdgePlots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

stats_file <- file.path(data_dir, "significant_edge_stats.csv")
stats_table <- read_csv(stats_file)

edge_files <- list.files(data_dir, pattern = "^edge_\\d+_to_\\d+\\.csv$", full.names = TRUE)

for (file in edge_files) {
  df <- read_csv(file)
  
  edge_name <- str_extract(basename(file), "\\d+_to_\\d+")
  from_to <- str_split(edge_name, "_to_")[[1]]
  from <- as.integer(from_to[1])
  to <- as.integer(from_to[2])
  
  row_stat <- stats_table[stats_table$from == from & stats_table$to == to, ]
  r_val <- round(row_stat$r, 3)
  p_val <- signif(row_stat$p, 3)
  
  df$edge_class <- ifelse(df$edge_strength == 0, "Zero", "Non-zero")
  
  p <- ggplot(df, aes(x = score, y = edge_strength)) +
    geom_point(
      aes(color = edge_class), 
      size = 1.2, alpha = 0.6, shape = 16
    ) +
    scale_color_manual(
      values = c("Zero" = "#D3D3D3", "Non-zero" = "#0072B2")
    ) +
    geom_smooth(
      data = df[df$edge_strength != 0, ],
      method = "lm", se = TRUE,
      color = "#222222", fill = "#B0B0B0",
      alpha = 0.25, linewidth = 0.5
    ) +
    labs(
      title = sprintf("Edge %d → %d", from, to),
      subtitle = sprintf("Spearman r = %.3f, p = %.3g", r_val, p_val),
      x = "Score",
      y = "Edge Strength"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
      axis.title = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
  
  filename_base <- paste0("edge_", edge_name)
  
  ggsave(file.path(plot_dir, paste0(filename_base, ".png")), p, width = 6, height = 5, dpi = 600)
  ggsave(file.path(plot_dir, paste0(filename_base, ".svg")), p, width = 6, height = 5)
  
  cat("🎯 已保存高质量图：", filename_base, "\n")
}