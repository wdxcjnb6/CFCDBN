# 加载必要包
library(R.matlab)
library(circlize)
library(RColorBrewer)
library(tools)
library(scales)

# 路径配置
base_path <- "E:/PAC_network/Result/Direct_net_result/Direct_net_figure/Data"
output_dir <- "E:/PAC_network/R_code/规模网络/Figure/"
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)

# .mat 文件列表
files <- c("Anxiety_group.mat", "HC_group.mat")
full_paths <- file.path(base_path, files)

# 定义功能区（顺序与连接矩阵一致）
regions <- c(rep("DMN", 10),
             rep("SN", 8),
             rep("FPN", 6),
             rep("DAN", 4),
             rep("SMN", 6))
region_names <- unique(regions)

# 可自定义中文或英文标签
region_labels <- c("默认网络DMN", "显著性网络SN", "前额叶网络FPN", "注意网络DAN", "运动网络SMN")
names(region_labels) <- region_names

# 颜色配置
network_colors <- c(
  "DMN" = "#424186B3",  # 深蓝紫（70%透明度）
  "SN"  = "#2B758EB3",  # 蓝绿
  "FPN" = "#20A486B3",  # 淡蓝灰
  "DAN" = "#6ECE58B3",  # 淡紫灰
  "SMN" = "#EFE51CB3"   # 柔和灰蓝紫
)

# 可自定义函数，用于绘制有向弦图
plot_directed_chord <- function(matrix_data, output_file) {
  n_regions <- length(region_names)
  agg_matrix <- matrix(0, nrow = n_regions, ncol = n_regions)
  rownames(agg_matrix) <- colnames(agg_matrix) <- region_names
  
  for (i in 1:n_regions) {
    for (j in 1:n_regions) {
      from_idx <- which(regions == region_names[i])
      to_idx <- which(regions == region_names[j])
      sub_matrix <- matrix_data[from_idx, to_idx]
      agg_matrix[i, j] <- mean(sub_matrix, na.rm = TRUE)
    }
  }
  
  # 转换为边数据框
  edge_df <- as.data.frame(as.table(agg_matrix))
  colnames(edge_df) <- c("from", "to", "value")
  edge_df <- edge_df[edge_df$value > 0, ]
  
  # 透明度映射
  edge_df$transparency <- rescale(edge_df$value, to = c(0.8, 0.3))
  
  # 颜色配置
  grid_col <- setNames(network_colors[region_names], region_names)
  
  # 封装绘图逻辑为函数，避免重复
  draw_chord <- function() {
    circos.clear()
    circos.par(message = FALSE, gap.after = c(rep(10, n_regions - 1), 8))
    
    chordDiagram(
      edge_df[, c("from", "to", "value")],
      grid.col = grid_col,
      directional = -1,
      direction.type = c("diffHeight", "arrows"),
      diffHeight = 0.05,
      transparency = edge_df$transparency,
      link.arr.type = "big.arrow",
      annotationTrack = "grid",
      link.sort = TRUE,
      link.decreasing = FALSE,
      link.lwd = 2,
      link.border = NA,
      preAllocateTracks = 1
    )
    
    circos.trackPlotRegion(
      track.index = 1,
      panel.fun = function(x, y) {
        sector_name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(
          mean(xlim),
          ylim[1] + 0.2,
          sector_name,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 1.3,
          font = 2
        )
      },
      bg.border = NA
    )
  }
  
  # ==== SVG 保存 ====
  svg_file <- paste0(output_file, ".svg")
  grDevices::svg(file = svg_file, width = 7, height = 7, onefile = TRUE, bg = "transparent")
  draw_chord()
  dev.off()
  
  # ==== PNG 保存 ====
  png_file <- paste0(output_file, ".png")
  png(
    filename = png_file,
    width = 1800,
    height = 1800,
    res = 300
  )
  draw_chord()
  dev.off()
}

# 循环读取.mat文件，生成相应的弦图
for (file in full_paths) {
  mat_data <- readMat(file)
  var_name <- names(mat_data)[1]
  conn_matrix <- mat_data[[var_name]]
  
  if (length(dim(conn_matrix)) == 3) {
    conn_matrix <- conn_matrix[, , 1]
  }
  
  output_file <- file.path(output_dir, paste0(file_path_sans_ext(basename(file)), "_chord"))
  
  plot_directed_chord(conn_matrix, output_file)
}
