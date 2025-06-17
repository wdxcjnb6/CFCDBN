library(ggplot2)
library(ggrain)
library(svglite)  # 确保安装，用于保存 SVG 图像

# 读取数据
hc_data <- read.csv("E:/PAC_network/Result/PAC_result/Fig/subject_pac_values_HC.csv")
anxiety_data <- read.csv("E:/PAC_network/Result/PAC_result/Fig/subject_pac_values_Anxiety.csv")

# 合并数据
df <- rbind(hc_data, anxiety_data)

# 设置配色
color_scheme <- c("HC" = "#1f77b4", "Anxiety" = "#d62728")

# 绘图
p <- ggplot(df, aes(x = ROIs, y = PAC_Value, fill = Group, color = Group)) +
  geom_rain(alpha = 0.5,
            boxplot.args = list(color = "black", outlier.shape = NA, size = 0.5)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  labs(x = "ROI", y = "PAC Value", title = "Distribution of HC and Anxiety PAC Values by ROI") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 保存路径
save_path <- "E:/PAC_network/R_code/Rain/Figure"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# 保存为 SVG
ggsave(filename = file.path(save_path, "PAC_group_distribution.svg"),
       plot = p, device = svglite::svglite, width = 10, height = 6)

cat("✅ SVG 图已成功保存至：", file.path(save_path, "PAC_group_distribution.svg"), "\n")
p