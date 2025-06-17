# 加载必要包
library(R.matlab)
library(circlize)
library(RColorBrewer)
library(tools)
library(scales)

# 路径配置
base_path <- "E:/PAC_network/Result/Direct_net_result/Direct_net_figure/Data"
output_dir <- "E:/PAC_network/R_code/规模网络/Figure/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 文件列表
files <- c("significant_links_group.mat")
full_paths <- file.path(base_path, files)

# 你提供的 34 个 ROI 名称
roi_names <- c(
  'HIP.L','HIP.R','SFGmed.L','SFGmed.R','PCUN.L','PCUN.R','ANG.L','ANG.R',
  'PCG.L','PCG.R','INS.L','INS.R','ACG.L','ACG.R','AMYG.L','AMYG.R',
  'THA.L','THA.R','MFG.L','MFG.R','IPL.L','IPL.R','IFGtriang.L','IFGtriang.R',
  'SFGdor.L','SFGdor.R','SPG.L','SPG.R','PoCG.L','PoCG.R','PreCG.L','PreCG.R',
  'SMA.L','SMA.R'
)

# 映射每个 ROI 所属的功能区（与聚合版颜色保持一致）
roi_to_network <- c(
  'HIP.L' = "DMN",  'HIP.R' = "DMN",
  'SFGmed.L' = "DMN",  'SFGmed.R' = "DMN",
  'PCUN.L' = "DMN",  'PCUN.R' = "DMN",
  'ANG.L' = "DMN",   'ANG.R' = "DMN",
  'PCG.L' = "DMN",   'PCG.R' = "DMN",
  
  'INS.L' = "SN",    'INS.R' = "SN",
  'ACG.L' = "SN",    'ACG.R' = "SN",
  'AMYG.L' = "SN",   'AMYG.R' = "SN",
  'THA.L' = "SN",    'THA.R' = "SN",
  
  'MFG.L' = "FPN",   'MFG.R' = "FPN",
  'IPL.L' = "FPN",   'IPL.R' = "FPN",
  'IFGtriang.L' = "FPN", 'IFGtriang.R' = "FPN",
  
  'SFGdor.L' = "DAN", 'SFGdor.R' = "DAN",
  'SPG.L' = "DAN",    'SPG.R' = "DAN",
  
  'PoCG.L' = "SMN",   'PoCG.R' = "SMN",
  'PreCG.L' = "SMN",  'PreCG.R' = "SMN",
  'SMA.L' = "SMN",    'SMA.R' = "SMN"
)

# 功能网络颜色（与你之前聚合网络图一致）
# 替换为你指定色条中提取的五种代表色
network_colors <- setNames(
  c(
    "DMN" = "#424186B3",  # 深蓝紫（70%透明度）
    "SN"  = "#2B758EB3",  # 蓝绿
    "FPN" = "#20A486B3",  # 淡蓝灰
    "DAN" = "#6ECE58B3",  # 淡紫灰
    "SMN" = "#EFE51CB3"   # 柔和灰蓝紫
  ),
  c("DMN", "SN", "FPN", "DAN", "SMN")
)

# 每个 ROI 分配颜色
roi_colors <- setNames(
  scales::alpha(network_colors[roi_to_network], 0.8),
  names(roi_to_network)
)



#network_colors <- setNames(brewer.pal(5, "Set1"), c("DMN", "SN", "FPN", "DAN", "SMN"))

# 每个 ROI 分配颜色
#roi_colors <- setNames( scales::alpha(network_colors[roi_to_network], 0.8),names(roi_to_network))  



# 绘图函数（不做聚合，单脑区显示）
plot_directed_chord_roi <- function(matrix_data, output_file) {
  rownames(matrix_data) <- colnames(matrix_data) <- roi_names
  
  # 转换为边数据框
  edge_df <- as.data.frame(as.table(matrix_data))
  colnames(edge_df) <- c("from", "to", "value")
  edge_df <- edge_df[edge_df$value > 0, ]
  
  
  # 映射透明度
  edge_df$transparency <- rescale(edge_df$value, to = c(0.8, 0.3))
  
  draw_chord <- function() {
    circos.clear()
    
    # 根据实际扇形数量设置 gap.after，避免报错
    used_rois <- unique(c(edge_df$from, edge_df$to))
    gap_vals <- rep(14, length(used_rois))  # 每个扇区后间隔12度
    circos.par(gap.after = gap_vals)
    
    chordDiagram(
      edge_df[, c("from", "to", "value")],
      grid.col = roi_colors,
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
    
    # 标签
    circos.trackPlotRegion(
      track.index = 1,
      panel.fun = function(x, y) {
        sector_name <- get.cell.meta.data("sector.index")
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        circos.text(
          mean(xlim), ylim[1] + 0.2,
          sector_name,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 0.8,
          font = 2
        )
      },
      bg.border = NA
    )
  }
  
  # 保存 EPS
  # ==== SVG 保存 ====
  svg_file <- paste0(output_file, ".svg")
  grDevices::svg(file = svg_file, width = 7, height = 7, onefile = TRUE, bg = "transparent")
  draw_chord()
  dev.off()
  
  
  # 保存 PNG
  png_file <- paste0(output_file, ".png")
  png(filename = png_file, width = 1800, height = 1800, res = 300)
  draw_chord()
  dev.off()

}

# 主执行
for (file in full_paths) {
  mat_data <- readMat(file)
  var_name <- names(mat_data)[1]
  conn_matrix <- mat_data[[var_name]]
  
  if (length(dim(conn_matrix)) == 3) {
    conn_matrix <- conn_matrix[,,1]
  }
  
  output_file <- file.path(output_dir, paste0(file_path_sans_ext(basename(file)), "_ROI_chord"))
  plot_directed_chord_roi(conn_matrix, output_file)
}
