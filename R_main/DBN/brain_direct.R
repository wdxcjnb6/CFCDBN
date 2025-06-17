library(brainconn)

# 加载 brainconn 自带 AAL90 atlas
data("aal90", package = "brainconn")

# Step 1：确保 ROI.Name 是字符型
aal90$ROI.Name <- as.character(aal90$ROI.Name)

# Step 2：你自定义的边列表（有向 + 权重）
edges <- data.frame(
  source = c('PreCG.R', 'PreCG.R', 'ANG.R', 'PCG.L', 'SPG.L',
             'INS.R', 'PCUN.L', 'SFGdor.L', 'ANG.L', 'PCG.L'),
  target = c('HIP.L', 'SFGmed.L', 'PCUN.R', 'PCUN.R', 'PCUN.R',
             'PCG.L', 'INS.R', 'THA.R', 'MFG.R', 'PreCG.L'),
  weight = c(0.3519, 0.3618, 0.3179, -0.3383, 0.3320,
             0.3536, -0.3184, -0.3345, -0.4045, 0.3218)
)

# Step 3：获取涉及到的所有节点名称，并清理空格
all_nodes <- trimws(unique(c(edges$source, edges$target)))

# Step 4：匹配这些脑区在 AAL90 表中的索引
node_rows <- match(all_nodes, aal90$ROI.Name)

# 检查匹配是否成功
if (any(is.na(node_rows))) {
  cat("❌ 以下脑区名未匹配成功：\n")
  print(all_nodes[is.na(node_rows)])
  stop("匹配失败，请检查脑区名称拼写")
} else {
  cat("✅ 所有脑区匹配成功！\n")
}

# Step 5：构建连接矩阵（子图）
n <- length(all_nodes)
conmat <- matrix(0, nrow = n, ncol = n)
rownames(conmat) <- colnames(conmat) <- all_nodes

# 将边表信息写入连接矩阵
for (i in 1:nrow(edges)) {
  src <- match(edges$source[i], all_nodes)
  tgt <- match(edges$target[i], all_nodes)
  conmat[src, tgt] <- edges$weight[i]
}

# Step 6：提取对应坐标和标签
coords <- aal90[node_rows, c("x.mni", "y.mni", "z.mni")]
labels <- aal90$ROI.Name[node_rows]

# Step 7：使用 brainconn 自定义绘图
library(brainconn)

# 加载 brainconn 自带 AAL90 atlas
data("aal90", package = "brainconn")

# Step 1：确保 ROI.Name 是字符型
aal90$ROI.Name <- as.character(aal90$ROI.Name)

# Step 2：你自定义的边列表（有向 + 权重）
edges <- data.frame(
  source = c('PreCG.R', 'PreCG.R', 'ANG.R', 'PCG.L', 'SPG.L',
             'INS.R', 'PCUN.L', 'SFGdor.L', 'ANG.L', 'PCG.L'),
  target = c('HIP.L', 'SFGmed.L', 'PCUN.R', 'PCUN.R', 'PCUN.R',
             'PCG.L', 'INS.R', 'THA.R', 'MFG.R', 'PreCG.L'),
  weight = c(0.3519, 0.3618, 0.3179, -0.3383, 0.3320,
             0.3536, -0.3184, -0.3345, -0.4045, 0.3218)
)

# Step 3：获取涉及到的所有节点名称，并清理空格
all_nodes <- trimws(unique(c(edges$source, edges$target)))

# Step 4：匹配这些脑区在 AAL90 表中的索引
node_rows <- match(all_nodes, aal90$ROI.Name)

# 检查匹配是否成功
if (any(is.na(node_rows))) {
  cat("❌ 以下脑区名未匹配成功：\n")
  print(all_nodes[is.na(node_rows)])
  stop("匹配失败，请检查脑区名称拼写")
} else {
  cat("✅ 所有脑区匹配成功！\n")
}

# Step 5：构建连接矩阵（子图）
n <- length(all_nodes)
conmat <- matrix(0, nrow = n, ncol = n)
rownames(conmat) <- colnames(conmat) <- all_nodes

# 将边表信息写入连接矩阵
for (i in 1:nrow(edges)) {
  src <- match(edges$source[i], all_nodes)
  tgt <- match(edges$target[i], all_nodes)
  conmat[src, tgt] <- edges$weight[i]
}

# Step 6：提取对应坐标和标签
coords <- aal90[node_rows, c("x.mni", "y.mni", "z.mni")]
labels <- aal90$ROI.Name[node_rows]

# Step 7：使用 brainconn 自定义绘图
# Step 7（正确方式）：构造自定义 atlas 并绘图

# 构造自定义 atlas 数据框
custom_atlas <- data.frame(
  ROI.Name = labels,
  x.mni = coords[, 1],
  y.mni = coords[, 2],
  z.mni = coords[, 3],
  stringsAsFactors = FALSE
)

# 把它作为名为 "custom" 的 atlas 放入全局环境
assign("custom", custom_atlas, envir = .GlobalEnv)

brainconn(
  atlas = "custom",
  conmat = conmat,
  node.size = 6,
  edge.alpha = 0.9,
  node.color = "skyblue",  # ✅ 全局颜色，推荐测试时使用
  view = "right"
)



