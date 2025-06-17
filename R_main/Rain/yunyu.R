library(ggplot2)
library(ggrain)

# 示例数据
set.seed(123)
df <- data.frame(
  Group = rep(c("Group1", "Group2", "Group3", "Group4", "Group5"), each = 100),
  Condition = rep(c("HC", "Anxiety"), times = 250),
  Value = c(
    rnorm(50, mean = 3.5, sd = 0.5),  # Group1 HC
    rnorm(50, mean = 4.0, sd = 0.5),  # Group1 Anxiety
    rnorm(50, mean = 3.7, sd = 0.5),  # Group2 HC
    rnorm(50, mean = 4.2, sd = 0.5),  # Group2 Anxiety
    rnorm(50, mean = 3.3, sd = 0.5),  # Group3 HC
    rnorm(50, mean = 4.1, sd = 0.5),  # Group3 Anxiety
    rnorm(50, mean = 3.8, sd = 0.5),  # Group4 HC
    rnorm(50, mean = 4.3, sd = 0.5),  # Group4 Anxiety
    rnorm(50, mean = 3.6, sd = 0.5),  # Group5 HC
    rnorm(50, mean = 4.4, sd = 0.5)   # Group5 Anxiety
  )
)

# 配色方案：HC 为绿色，Anxiety 为蓝色
color_scheme <- c("HC" = "#2ca02c", "Anxiety" = "#1f77b4")  # 互换颜色

# 绘制分布图
ggplot(df, aes(x = Group, y = Value, fill = Condition, color = Condition)) +
  geom_rain(alpha = 0.5,  # 设置透明度为 0.5
            boxplot.args = list(color = "black", outlier.shape = NA, size = 0.5)) +  # 统一箱线图的颜色和边框
  theme_classic() +
  scale_fill_manual(values = color_scheme) +  # 使用统一的颜色方案
  scale_color_manual(values = color_scheme) +  # 统一散点图和箱线图的颜色
  labs(x = "Group", y = "Value", title = "Distribution of HC and Anxiety in Five Groups")
