# 加载必要包
library(circlize)
library(scales)  # 加载scales包，提供rescale函数

# 清空环形布局设置（关键步骤！）
circos.clear()

# 创建连接矩阵（这里只保留 t_stat 大于 1 的连接）
chord_matrix <- matrix(c(
  1.539045, 0.643067, 0.849735, 0.536657, 2.438411,  # DMN
  0.498035, 1.024674, 0.449856, 2.468088, 0.615169,  # SN
  1.070298, 0.314322, 0.912501, 0.341047, 0.160667,  # FPN
  0.121904, 0.279232, 0.417424, 1.085960, 0.812428,  # DAN
  0.086809, 0.871687, 2.842483, 0.909732, 0.329435   # SMN
), nrow = 5, ncol = 5, byrow = TRUE)

# 设置阈值，保留 t_stat > 2 的连接
chord_matrix[chord_matrix <= 2] <- 1 # 小于等于2的连接置为0

# 定义网络颜色（带透明度）
network_colors <- c(
  "DMN" = "#424186B3",  # 深蓝紫（70%透明度）
  "SN"  = "#2B758EB3",  # 蓝绿
  "FPN" = "#20A486B3",  # 淡蓝灰
  "DAN" = "#6ECE58B3",  # 淡紫灰
  "SMN" = "#EFE51CB3"   # 柔和灰蓝紫
)

# 设置行和列名为网络名称
rownames(chord_matrix) <- colnames(chord_matrix) <- names(network_colors)

# 创建一个边数据框并根据值进行透明度映射
edge_df <- as.data.frame(as.table(chord_matrix))
colnames(edge_df) <- c("from", "to", "value")
edge_df <- edge_df[edge_df$value > 0, ]

# 透明度映射
edge_df$transparency <- rescale(edge_df$value, to = c(0.8, 0.3))

# 颜色配置
grid_col <- setNames(network_colors[names(network_colors)], names(network_colors))

# 封装绘图逻辑为函数，避免重复
draw_chord <- function() {
  circos.clear()
  
  # 设置环形布局参数
  circos.par(message = FALSE, gap.after = c(rep(10, 4), 15))
  
  chordDiagram(
    edge_df[, c("from", "to", "value")],
    grid.col = grid_col,
    directional = -1,
    direction.type = c("diffHeight", "arrows"),
    diffHeight = 0.05,
    transparency = edge_df$transparency,
    link.arr.type = "big.arrow",  # 使用大箭头
    annotationTrack = "grid",
    link.sort = TRUE,
    link.decreasing = FALSE,
    link.lwd = 2,
    link.border = NA,
    preAllocateTracks = 1
  )
  
  # 添加旋转标签（优化显示）
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector_name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(mean(xlim), ylim[1] + 0.2,
                  sector_name,
                  facing = "clockwise",
                  niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = 1.3,
                  font = 2)
    },
    bg.border = NA
  )
}

# 保存图像
png("E:/PAC_network/R_code/规模网络/Figure2/Significant_Connections_Chord_Diagram.png", 
    width = 1200, height = 1000, res = 150)
svg("E:/PAC_network/R_code/规模网络/Figure2/Significant_Connections_Chord_Diagram.svg", 
    width = 7, height = 7, onefile = TRUE, bg = "transparent")

draw_chord()

dev.off()
circos.clear()

cat("🎯 已生成带方向指示的大箭头弦图：Significant_Connections_Chord_Diagram.png\n")
