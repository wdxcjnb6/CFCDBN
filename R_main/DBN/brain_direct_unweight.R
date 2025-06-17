library(brainconn)
library(ggplot2)
library(ggraph)
library(igraph)

# 加载 AAL90 atlas
data("aal90", package = "brainconn")
aal90$ROI.Name <- as.character(aal90$ROI.Name)

# ------------------【1. 读取 CSV 文件】------------------
# 读取 net1 和 net2 文件
net1_file <- "E:/PAC_network/R_code/有向脑网络/net1.csv"
net2_file <- "E:/PAC_network/R_code/有向脑网络/net2.csv"

# 读取文件
net1_edges <- read.csv(net1_file, stringsAsFactors = FALSE)
net2_edges <- read.csv(net2_file, stringsAsFactors = FALSE)

# ------------------【2. 确认列名是否为 Source 和 Target】------------------
colnames(net1_edges) <- c("source", "target")
colnames(net2_edges) <- c("source", "target")

# ------------------【3. 查看前几条边】------------------
head(net1_edges)
head(net2_edges)

# 函数：绘制脑网络图
plot_brain_network <- function(edges, plot_dir, atlas_name) {
  all_nodes <- unique(c(edges$source, edges$target))
  node_rows <- match(all_nodes, aal90$ROI.Name)
  coords_all <- aal90[node_rows, c("x.mni", "y.mni", "z.mni"), drop = FALSE]
  labels_all <- aal90$ROI.Name[node_rows]
  
  # 构造连接矩阵
  n_all <- length(all_nodes)
  conmat_all <- matrix(0, nrow = n_all, ncol = n_all)
  rownames(conmat_all) <- colnames(conmat_all) <- all_nodes
  for (i in 1:nrow(edges)) {
    src <- match(edges$source[i], all_nodes)
    tgt <- match(edges$target[i], all_nodes)
    conmat_all[src, tgt] <- 1
  }
  
  custom_atlas_all <- data.frame(
    ROI.Name = labels_all,
    x.mni = coords_all[, 1],
    y.mni = coords_all[, 2],
    z.mni = coords_all[, 3],
    stringsAsFactors = FALSE
  )
  assign("custom", custom_atlas_all, envir = .GlobalEnv)
  
  # 定义视角
  views <- c("top", "bottom", "left", "right", "front", "back")
  
  # 循环绘制不同视角
  for (view in views) {
    p <- brainconn(
      atlas = "custom",
      conmat = conmat_all,
      node.size = 5,
      edge.width = 1,
      edge.alpha = 0.8,
      node.color = "skyblue",
      view = view
    ) + geom_node_text(aes(label = name), size = 3, repel = TRUE,
                       fontface = "bold", color = "black")
    
    # 保存图形
    save_dir <- file.path(plot_dir, atlas_name)
    if (!dir.exists(save_dir)) dir.create(save_dir)
    
    # 保存为 PNG 和 EPS，增加图像质量（dpi = 600）
    ggsave(file.path(save_dir, paste0("p_", view, ".png")), p, width = 8, height = 6, dpi = 600)
    ggsave(file.path(save_dir, paste0("p_", view, ".svg")), p, device = cairo_ps, width = 8, height = 6)
    
    cat("🎯 已保存图：", atlas_name, "视角：", view, "\n")
  }
}

# ====================
# 绘制 net1 和 net2 的图
# ====================
plot_dir <- "E:/PAC_network/R_code/有向脑网络/Figure"

# 绘制 net1 图
plot_brain_network(net1_edges, plot_dir, "net1")

# 绘制 net2 图
plot_brain_network(net2_edges, plot_dir, "net2")
