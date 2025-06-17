# ========== 自动安装并加载必要包 ==========
required_packages <- c("ggplot2", "readr", "svglite")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg)
  }
}

library(ggplot2)
library(readr)

# ========== Set Paths ==========
data_dir <- "E:/PAC_network/R_code/PAC_strength/Data"
fig_dir  <- "E:/PAC_network/R_code/PAC_strength/Figure"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ========== Define Colors ==========
group_colors <- c("Anxiety" = "#27ae60", "HC" = "#27ae60")
network_palette <- c(
  "DMN" = "#3E4F94",
  "SN"  = "#3E90BF",
  "FPN" = "#A6C0E3",
  "DAN" = "#BCA9D9",
  "SMN" = "#A8B0C9"
)

# ========== ROI → Network Mapping ==========
roi_order <- c(
  "HIP.L", "HIP.R", "SFGmed.L", "SFGmed.R", "PCUN.L", "PCUN.R", "ANG.L", "ANG.R", "PCG.L", "PCG.R",
  "INS.L", "INS.R", "ACG.L", "ACG.R", "AMYG.L", "AMYG.R", "THA.L", "THA.R",
  "MFG.L", "MFG.R", "IPL.L", "IPL.R", "IFGtriang.L", "IFGtriang.R",
  "SFGdor.L", "SFGdor.R", "SPG.L", "SPG.R",
  "PoCG.L", "PoCG.R", "PreCG.L", "PreCG.R", "SMA.L", "SMA.R"
)
roi_network_map <- data.frame(
  ROI = roi_order,
  Network = c(rep("DMN", 10),
              rep("SN", 8),
              rep("FPN", 6),
              rep("DAN", 4),
              rep("SMN", 6))
)

# ========== Save Plot (SVG + PNG) ==========
save_plot_safe <- function(plot, filename_base) {
  svg_file <- file.path(fig_dir, paste0(filename_base, ".svg"))
  png_file <- file.path(fig_dir, paste0(filename_base, ".png"))
  
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(svg_file, plot = plot, width = 12, height = 5, device = "svg", dpi = 600)
    message("✅ SVG saved: ", svg_file)
  } else {
    warning("⚠️ svglite not installed, skipping SVG.")
  }
  ggsave(png_file, plot = plot, width = 12, height = 5, dpi = 600)
  message("✅ PNG saved: ", png_file)
}

# ========== Main Loop for Both Groups ==========
for (group in c("Anxiety", "HC")) {
  
  # === Read Data ===
  df <- read_csv(file.path(data_dir, paste0(group, "_PAC_sorted_for_R.csv")))
  
  # === Merge with Network Info ===
  df <- merge(df, roi_network_map, by = "ROI", sort = FALSE)
  df$ROI <- factor(df$ROI, levels = roi_order)  # fixed x-axis order
  df$Index <- as.integer(factor(df$ROI, levels = roi_order))
  
  # === Shrink Error Bars ===
  n <- 40
  se_scale <- 0.5
  df$SE <- df$Std / sqrt(n)
  df$ymin <- df$Mean - df$SE * se_scale
  df$ymax <- df$Mean + df$SE * se_scale
  
  # === Plot ===
  p <- ggplot(df, aes(x = Index, y = Mean)) +
    geom_line(color = group_colors[group], linewidth = 0.9) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax, color = Network),
                  width = 0.15, linewidth = 0.6) +
    geom_point(aes(color = Network), size = 3.2) +
    scale_x_continuous(breaks = df$Index, labels = df$ROI) +
    scale_color_manual(
      values = network_palette,
      breaks = names(network_palette)  # 👈 强制图例顺序匹配定义顺序
    ) +
    labs(
      title = paste0(group, " group: PAC trend by ROI"),
      x = "Brain region",
      y = "Mean PAC",
      color = "Functional network"
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major.y = element_line(color = "gray85", linetype = "dashed"),
      legend.position = "top"
    )
  
  # === Save ===
  save_plot_safe(p, paste0(group, "_PAC_by_network"))
}
p