# ------------------【加载包】------------------
library(ggplot2)
library(ggrain)      # 安装：remotes::install_github("kassambara/ggrain")
library(readr)

# ------------------【读取数据】------------------
df <- read_csv("E:/PAC_network/R_code/fooof/beta_individual_peak_data.csv")

# 确保 Group 为因子并清除空格
df$Group <- trimws(df$Group)
df$Group <- factor(df$Group, levels = c("HC", "Anxiety"))  # 顺序：HC 在前，Anxiety 在后

# ------------------【定义颜色方案】------------------
color_scheme <- c("HC" = "#2ca02c",       # 绿色
                  "Anxiety" = "#1f77b4")  # 蓝色

# ------------------【保存路径】------------------
save_path <- "E:/PAC_network/R_code/fooof/Figure"
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# ------------------【1. 按组分开数据并计算均值和方差】------------------
hc_power <- df$Beta_Peak_Power[df$Group == "HC"]
anxiety_power <- df$Beta_Peak_Power[df$Group == "Anxiety"]

hc_freq <- df$Beta_Peak_Frequency[df$Group == "HC"]
anxiety_freq <- df$Beta_Peak_Frequency[df$Group == "Anxiety"]

hc_bw <- df$Beta_BW[df$Group == "HC"] # 获取 HC 组的 Beta 带宽数据
anxiety_bw <- df$Beta_BW[df$Group == "Anxiety"] # 获取 Anxiety 组的 Beta 带宽数据

# 计算均值和方差
mean_hc_power <- mean(hc_power)
var_hc_power <- var(hc_power)

mean_anxiety_power <- mean(anxiety_power)
var_anxiety_power <- var(anxiety_power)

mean_hc_freq <- mean(hc_freq)
var_hc_freq <- var(hc_freq)

mean_anxiety_freq <- mean(anxiety_freq)
var_anxiety_freq <- var(anxiety_freq)

mean_hc_bw <- mean(hc_bw)
var_hc_bw <- var(hc_bw)

mean_anxiety_bw <- mean(anxiety_bw)
var_anxiety_bw <- var(anxiety_bw)

# ------------------【2. 独立样本 t 检验】------------------
# 峰值功率比较
t_power <- t.test(hc_power, anxiety_power)
# 峰值频率比较
t_freq <- t.test(hc_freq, anxiety_freq)
# 带宽比较
t_bw <- t.test(hc_bw, anxiety_bw)

# 输出 t 检验结果
cat("✅ β频段峰值功率：t = ", round(t_power$statistic, 3), ", p = ", round(t_power$p.value, 4), "\n")
cat("   HC 组均值 = ", round(mean_hc_power, 4), ", 方差 = ", round(var_hc_power, 4), "\n")
cat("   Anxiety 组均值 = ", round(mean_anxiety_power, 4), ", 方差 = ", round(var_anxiety_power, 4), "\n")

cat("✅ β频段峰值频率：t = ", round(t_freq$statistic, 3), ", p = ", round(t_freq$p.value, 4), "\n")
cat("   HC 组均值 = ", round(mean_hc_freq, 4), ", 方差 = ", round(var_hc_freq, 4), "\n")
cat("   Anxiety 组均值 = ", round(mean_anxiety_freq, 4), ", 方差 = ", round(var_anxiety_freq, 4), "\n")

cat("✅ β频段带宽：t = ", round(t_bw$statistic, 3), ", p = ", round(t_bw$p.value, 4), "\n")
cat("   HC 组均值 = ", round(mean_hc_bw, 4), ", 方差 = ", round(var_hc_bw, 4), "\n")
cat("   Anxiety 组均值 = ", round(mean_anxiety_bw, 4), ", 方差 = ", round(var_anxiety_bw, 4), "\n")

# ------------------【3. 绘制 Beta 峰值功率分布图】------------------
p1 <- ggplot(df, aes(x = Group, y = Beta_Peak_Power, fill = Group, color = Group)) +
  geom_rain(alpha = 0.5,
            boxplot.args = list(color = "black", outlier.shape = NA, size = 0.5)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  labs(x = "Group", y = "Beta Peak Power",
       title = "Distribution of Beta Peak Power by Group") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

# 保存 Beta 峰值功率图
ggsave(file.path(save_path, "Beta_Peak_Power_Distribution.eps"), plot = p1, device = cairo_ps, width = 8, height = 6)
ggsave(file.path(save_path, "Beta_Peak_Power_Distribution.png"), plot = p1, device = "png", dpi = 300, width = 8, height = 6)
ggsave(file.path(save_path, "Beta_Peak_Power_Distribution.svg"), plot = p1, device = "svg", width = 8, height = 6)

# ------------------【4. 绘制 Beta 峰值频率分布图】------------------
p2 <- ggplot(df, aes(x = Group, y = Beta_Peak_Frequency, fill = Group, color = Group)) +
  geom_rain(alpha = 0.5,
            boxplot.args = list(color = "black", outlier.shape = NA, size = 0.5)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  labs(x = "Group", y = "Beta Peak Frequency (Hz)",
       title = "Distribution of Beta Peak Frequency by Group") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

# 保存 Beta 峰值频率图
ggsave(file.path(save_path, "Beta_Peak_Frequency_Distribution.eps"), plot = p2, device = cairo_ps, width = 8, height = 6)
ggsave(file.path(save_path, "Beta_Peak_Frequency_Distribution.png"), plot = p2, device = "png", dpi = 300, width = 8, height = 6)
ggsave(file.path(save_path, "Beta_Peak_Frequency_Distribution.svg"), plot = p2, device = "svg", width = 8, height = 6)

# ------------------【5. 绘制 Beta 带宽分布图】------------------
p3 <- ggplot(df, aes(x = Group, y = Beta_BW, fill = Group, color = Group)) +
  geom_rain(alpha = 0.5,
            boxplot.args = list(color = "black", outlier.shape = NA, size = 0.5)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  labs(x = "Group", y = "Beta Bandwidth (Hz)",
       title = "Distribution of Beta Bandwidth by Group") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

# 保存 Beta 带宽图
ggsave(file.path(save_path, "Beta_Bandwidth_Distribution.eps"), plot = p3, device = cairo_ps, width = 8, height = 6)
ggsave(file.path(save_path, "Beta_Bandwidth_Distribution.png"), plot = p3, device = "png", dpi = 300, width = 8, height = 6)
ggsave(file.path(save_path, "Beta_Bandwidth_Distribution.svg"), plot = p3, device = "svg", width = 8, height = 6)

