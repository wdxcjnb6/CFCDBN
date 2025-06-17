library(ggplot2)
library(ggdist)
library(dplyr)

# 读取数据
hc_data <- read.csv("E:/PAC_network/Result/PAC_result/Fig/subject_pac_values_HC.csv")
anxiety_data <- read.csv("E:/PAC_network/Result/PAC_result/Fig/subject_pac_values_Anxiety.csv")

# 给数据添加分组标签
hc_data$Group <- "HC"
anxiety_data$Group <- "Anxiety"

# 合并数据
df <- rbind(hc_data, anxiety_data)

# 设置配色
color_scheme <- c("HC" = "#AB97E0", "Anxiety" = "#F58F96")

# 绘图
p <- ggplot(df, aes(x = ROIs, y = PAC_Value, fill = Group)) +
  
  # 半眼云图：Group 分组，重叠显示
  ggdist::stat_halfeye(
    aes(fill = Group),
    adjust = 0.5,
    justification = -.3,
    .width = 0,
    point_colour = NA,
    alpha = 0.9
  ) +
  
  # 箱线图：Group 分组，微小偏移
  geom_boxplot(
    aes(color = Group),
    width = 0.12,
    outlier.color = NA,
    alpha = 0.5,
    position = position_dodge(width = 0.2)  # 控制偏移距离
  ) +
  
  # 散点图：不偏移，跟云图对齐
  ggdist::stat_dots(
    aes(color = Group),
    side = "left",
    justification = 1.1,
    binwidth = 0.0015,  # 更小的 binwidth
    size = 0.8,  # 减小散点大小
    alpha = 0.4,  # 增加透明度
    overflow = "compress"  # 压缩散点间距，避免溢出
  ) +
  
  # 颜色设置
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")  # 增加每行之间的间距
  ) +
  labs(
    title = "HC 和 Anxiety PAC 值在不同 ROI 中的分布",
    x = "ROI",
    y = "PAC 值"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存为 PNG 格式
ggsave("E:/PAC_network/R_code/Rain/Figure/PAC_distribution_HC_Anxiety.png", p, 
       width = 10, height = 8, dpi = 300)

# 保存为 SVG 格式
ggsave("E:/PAC_network/R_code/Rain/Figure/PAC_distribution_HC_Anxiety.svg", p, 
       width = 10, height = 8)
p