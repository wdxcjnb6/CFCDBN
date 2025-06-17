library(ggplot2)
library(patchwork)
library(reshape2)

# 读取数据
data <- read.csv("E:/PAC_network/R_code/leidatu/left_right_radar_data.csv")

# 将数据转换为长格式
data_long <- melt(data, id.vars = "network_names", 
                  measure.vars = c("left_mean_anxiety", "left_mean_hc", "right_mean_anxiety", "right_mean_hc"),
                  variable.name = "Group", value.name = "Mean")

# 处理数据，使左侧条形图数值取负
data_long$Mean[data_long$Group %in% c("left_mean_anxiety", "left_mean_hc")] <- 
  -data_long$Mean[data_long$Group %in% c("left_mean_anxiety", "left_mean_hc")]

# 添加标准差列
data_long$Std[data_long$Group == "left_mean_anxiety"] <- data$left_std_anxiety
data_long$Std[data_long$Group == "left_mean_hc"] <- data$left_std_hc
data_long$Std[data_long$Group == "right_mean_anxiety"] <- data$right_std_anxiety
data_long$Std[data_long$Group == "right_mean_hc"] <- data$right_std_hc

# 设置颜色
color_palette <- c("left_mean_anxiety" = "#FF6600", "left_mean_hc" = "#FF9966",
                   "right_mean_anxiety" = "#FF3300", "right_mean_hc" = "#FF6600")

# 绘制左侧条形图（左脑部分）
left_plot <- ggplot(data_long[data_long$Group %in% c("left_mean_anxiety", "left_mean_hc"), ], 
                    aes(x = network_names, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "white") + 
  coord_flip() + 
  labs(x = "Network", y = "Mean Value") +  # 添加横坐标标签
  scale_fill_manual(values = color_palette) +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2, position = position_dodge(width = 0.6)) +
  theme_minimal() +
  guides(fill = FALSE)  # 移除图例

# 绘制右侧条形图（右脑部分）
right_plot <- ggplot(data_long[data_long$Group %in% c("right_mean_anxiety", "right_mean_hc"), ], 
                     aes(x = network_names, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", widths = 0.6, color = "white") + 
  coord_flip() + 
  scale_fill_manual(values = color_palette) +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2, position = position_dodge(width = 0.6)) +
  scale_x_discrete(position = "top") +  # 设置右脑的横坐标位置
  labs(x = "Network", y = "Mean Value") +  # 添加横坐标标签
  theme_minimal() +
  guides(fill = FALSE)  # 移除图例

# 合并两个图表
final_plot <- left_plot + right_plot

# 保存图表为PNG格式
ggsave("E:/PAC_network/R_code/leidatu/final_plot.png", final_plot, width = 12, height = 8, dpi = 300)

# 保存图表为SVG格式
ggsave("E:/PAC_network/R_code/leidatu/final_plot.svg", final_plot, width = 12, height = 8)
