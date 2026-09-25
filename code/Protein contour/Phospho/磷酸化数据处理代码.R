library(tidyverse)
library(readxl)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# ==============================================================================
# 1. 数据准备与计算核心
# ==============================================================================
data_path <- "D:\\博士\\Phospho\\Full_data.xlsx" # 修改为你的实际路径

cat(">>> 正在读取数据...\n")
df_prot <- read_excel(data_path, sheet = "S5B-Log2 proc HeLa+EGF PROT")
df_phos <- read_excel(data_path, sheet = "S5D-Log2 proc HeLa+EGF PHOS")
colnames(df_prot)[1] <- "Gene"
colnames(df_phos)[1] <- "Phos_Key"

# --- 辅助函数：计算主要定位 (Main Fraction) ---
get_main_location <- function(df, time_pattern) {
  sub_df <- df %>% dplyr::select(contains(time_pattern))
  if(ncol(sub_df) != 24) return(NULL)
  
  mat_data <- 2^(as.matrix(sub_df)) # 反log
  
  # 计算每个Fraction的平均值
  # 结构: FR1, FR1, FR1, FR1, FR2...
  group_factor <- rep(paste0("FR", 1:6), each = 4)
  
  row_means <- t(apply(mat_data, 1, function(x) {
    tapply(x, group_factor, mean, na.rm=TRUE)
  }))
  
  # 找出最大值所在的 Fraction (即主要定位)
  # 并计算最大值的占比 (用于过滤低质量数据)
  results <- data.frame(
    ID = df[[1]],
    Main_Loc = apply(row_means, 1, function(x) {
      if(sum(x, na.rm=TRUE) == 0) return(NA)
      names(which.max(x)) # 返回 "FR1", "FR4" 等
    }),
    Max_Prop = apply(row_means, 1, function(x) {
      s <- sum(x, na.rm=TRUE)
      if(s == 0) return(0)
      max(x, na.rm=TRUE) / s
    }),
    stringsAsFactors = FALSE
  )
  return(results)
}

# --- 主循环：计算所有时间点的空间偏好 ---
timepoints <- c("EGF_CTRL", "EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min")
all_bias_data <- list()

cat(">>> 开始计算所有时间点的空间偏好...\n")

for (tp in timepoints) {
  # 1. 算总蛋白主要在哪里
  loc_prot <- get_main_location(df_prot, tp)
  if(is.null(loc_prot)) next
  
  # [修复步骤]: 总蛋白表的 ID 其实就是 Gene Name。
  # 我们把 loc_prot 的第一列 "ID" 改名为 "Gene"，这样 inner_join 就能找到了。
  colnames(loc_prot)[colnames(loc_prot) == "ID"] <- "Gene"
  
  # 2. 算磷酸化主要在哪里
  loc_phos <- get_main_location(df_phos, tp)
  if(is.null(loc_phos)) next
  
  # 3. 匹配
  # 从磷酸化 Key (如 AAK1_S623_M2) 提取基因名 (AAK1)
  loc_phos$Gene <- sub("_.*", "", loc_phos$ID)
  
  # 现在 loc_phos 和 loc_prot 都有 "Gene" 列了，可以合并
  merged <- inner_join(loc_phos, loc_prot, by = "Gene", suffix = c("_Phos", "_Prot")) %>%
    filter(!is.na(Main_Loc_Phos) & !is.na(Main_Loc_Prot)) %>%
    filter(Max_Prop_Phos > 0.4) # 过滤：只保留定位比较清晰的位点
  
  # 4. 定义"偏好" (Bias)
  bias_events <- merged %>%
    mutate(
      Timepoint = tp,
      Is_Biased = (Main_Loc_Phos != Main_Loc_Prot) # 定位不一致即为Bias
    )
  
  all_bias_data[[tp]] <- bias_events
}

# 合并大表
full_bias_df <- bind_rows(all_bias_data)
cat(">>> 计算完成！共分析了", nrow(full_bias_df), "个配对数据。\n")


# ==============================================================================
# 2. 可视化方案 A：网格热图 (Grid Plot with Numbers)
# ==============================================================================
# 目标：画一张图，展示某个时间点（如 20min）的转移矩阵

plot_bias_grid <- function(data_df, timepoint_name) {
  
  # 1. 筛选特定时间点的数据
  subset_df <- data_df %>% 
    filter(Timepoint == timepoint_name)
  
  # 2. 统计转换数量 (Prot Loc -> Phos Loc)
  # 这里的逻辑是：总蛋白在 Row，磷酸化在 Col
  trans_matrix <- table(subset_df$Main_Loc_Prot, subset_df$Main_Loc_Phos)
  
  # 补全可能缺失的 FR (确保矩阵是 6x6)
  all_fr <- paste0("FR", 1:6)
  mat_final <- matrix(0, nrow=6, ncol=6, dimnames = list(all_fr, all_fr))
  
  # 填入数据
  rows <- rownames(trans_matrix)
  cols <- colnames(trans_matrix)
  mat_final[rows, cols] <- trans_matrix
  
  # 3. 绘图
  # 对角线上的数字（Bias=0）通常很大，可能会掩盖非对角线的细节
  # 我们可以把对角线设为灰色或者保留，看你需求。这里保留但用不同颜色。
  
  col_fun = colorRamp2(c(0, max(mat_final)*0.3, max(mat_final)), c("white", "orange", "firebrick"))
  
  Heatmap(mat_final,
          name = "Count",
          col = col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          rect_gp = gpar(col = "white", lwd = 2), # 格子边框
          column_title = paste0("Phospho-site Location (", timepoint_name, ")"),
          row_title = "Total Protein Location",
          cell_fun = function(j, i, x, y, width, height, fill) {
            # 在格子里写数字
            grid.text(sprintf("%d", mat_final[i, j]), x, y, gp = gpar(fontsize = 10))
          })
}

# 测试：画一下 20min 的网格图
cat(">>> 正在绘制网格图 (Grid Plot)...\n")
draw(plot_bias_grid(full_bias_df, "EGF_20min"))


# ==============================================================================
# 3. 可视化方案 B：空间偏好网络图 (Spatial Bias Network)
# ==============================================================================
# 目标：用箭头展示“总蛋白 -> 磷酸化”的偏离趋势

plot_bias_network <- function(data_df, timepoint_name) {
  
  # 1. 准备数据
  # 我们只关心发生了"偏离"的那些 (去掉对角线)
  edges_df <- data_df %>%
    filter(Timepoint == timepoint_name) %>%
    filter(Main_Loc_Prot != Main_Loc_Phos) %>% # 只画不一致的
    group_by(Main_Loc_Prot, Main_Loc_Phos) %>%
    summarise(weight = n(), .groups = "drop") %>%
    filter(weight > 2) # 过滤掉只有1-2个位点的偶然事件，图会更干净
  
  if(nrow(edges_df) == 0) {
    message("没有检测到显著的偏离事件。")
    return(NULL)
  }
  
  # 2. 定义节点 (固定的细胞结构)
  nodes <- data.frame(
    id = paste0("FR", 1:6),
    label = c("FR1\nCytosol", "FR2\nRibosome", "FR3\nER/Golgi", 
              "FR4\nMitochondria", "FR5\nNuc.Peri", "FR6\nNucleus"),
    group = c("Cyto", "Cyto", "Memb", "Memb", "Nuc", "Nuc"),
    stringsAsFactors = FALSE
  )
  
  # 3. 构建 Graph
  g <- graph_from_data_frame(d = edges_df, vertices = nodes, directed = TRUE)
  
  # 4. 设置视觉属性
  
  # 颜色: 胞质绿，膜橙，核蓝
  node_cols <- c("Cyto"="#66C2A5", "Memb"="#FC8D62", "Nuc"="#8DA0CB")
  V(g)$color <- node_cols[V(g)$group]
  V(g)$frame.color <- "white"
  V(g)$label.color <- "black"
  V(g)$label.cex <- 0.9
  V(g)$size <- 35 # 节点大小
  
  # 边: 粗细代表数量
  E(g)$width <- sqrt(E(g)$weight) * 1.5 
  E(g)$arrow.size <- 0.8
  E(g)$color <- adjustcolor("gray40", alpha.f = 0.7)
  
  # 关键：给边加上数字标签 (你想要的数字！)
  E(g)$label <- E(g)$weight
  E(g)$label.cex <- 0.8
  E(g)$label.color <- "black"
  E(g)$label.dist <- 0 # 标签在正中间
  
  # 5. 布局 (手动设置坐标，模拟细胞结构)
  # 左上(FR1), 右上(FR2)
  # 左中(FR3), 右中(FR4)
  # 左下(FR5), 右下(FR6)
  layout_coords <- matrix(c(
    -1,  1,  # FR1
    1,  1,  # FR2
    -1,  0,  # FR3
    1,  0,  # FR4
    -1, -1,  # FR5
    1, -1   # FR6
  ), ncol = 2, byrow = TRUE)
  
  # 6. 绘图
  par(mar = c(1, 1, 3, 1))
  plot(g, 
       layout = layout_coords,
       edge.curved = 0.2, # 弧线，防止来回箭头重叠
       main = paste("Spatial Bias Network:", timepoint_name, "\n(Arrow: Protein Loc -> Phospho Loc)"))
  
  legend("topright", 
         legend = c("Cytoplasm", "Membrane", "Nucleus"), 
         fill = unique(node_cols), 
         bty = "n", cex = 0.8)
}

# 测试：画一下 20min 的网络图
cat(">>> 正在绘制网络图 (Network Plot)...\n")
plot_bias_network(full_bias_df, "EGF_20min")



library(tidyverse)
library(ggplot2)

# ==============================================================================
# 进阶方案：相对于 Control 的变化热图 (Delta Heatmap)
# ==============================================================================

plot_delta_grid <- function(bias_df) {
  
  # 1. 基础数据统计
  # 统计每个时间点、每个位置对的数量
  counts_data <- bias_df %>%
    group_by(Timepoint, Main_Loc_Prot, Main_Loc_Phos) %>%
    summarise(Count = n(), .groups = "drop")
  
  # 2. 补全数据 (FR1-FR6 全排列)
  all_fractions <- paste0("FR", 1:6)
  
  counts_complete <- counts_data %>%
    complete(Timepoint, 
             Main_Loc_Prot = all_fractions, 
             Main_Loc_Phos = all_fractions, 
             fill = list(Count = 0))
  
  # 3. 分离 Control 组作为基准
  ctrl_data <- counts_complete %>%
    filter(Timepoint == "EGF_CTRL") %>%
    select(Main_Loc_Prot, Main_Loc_Phos, Count_Ctrl = Count)
  
  # 4. 计算差异 (Delta)
  # 我们只保留处理组 (2min, 8min, 20min, 90min)
  plot_data <- counts_complete %>%
    filter(Timepoint != "EGF_CTRL") %>%
    left_join(ctrl_data, by = c("Main_Loc_Prot", "Main_Loc_Phos")) %>%
    mutate(
      # 计算差值：处理组 - 对照组
      Delta_Count = Count - Count_Ctrl,
      
      # 为了画图好看，我们可以根据差值生成标签
      # 只有变化量超过一定阈值（比如变动了2个以上）才显示数字，避免图太乱
      Label = ifelse(abs(Delta_Count) >= 2, sprintf("%+d", Delta_Count), "")
    ) %>%
    # 调整因子顺序
    mutate(
      Main_Loc_Prot = factor(Main_Loc_Prot, levels = rev(all_fractions)),
      Main_Loc_Phos = factor(Main_Loc_Phos, levels = all_fractions),
      Timepoint = factor(Timepoint, levels = c("EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min"))
    )
  
  # 5. 绘图 (红蓝配色：红增蓝减)
  # 设置颜色范围极限，保证0在中间是白色
  max_val <- max(abs(plot_data$Delta_Count), na.rm = TRUE)
  limit_range <- c(-max_val, max_val)
  
  p <- ggplot(plot_data, aes(x = Main_Loc_Phos, y = Main_Loc_Prot)) +
    # 画格子
    geom_tile(aes(fill = Delta_Count), color = "grey80", size = 0.3) +
    
    # 填入变化数值 (带正负号)
    geom_text(aes(label = Label), size = 3, color = "black") +
    
    # 分面展示
    facet_grid(. ~ Timepoint) +
    
    # 关键：红白蓝配色 (Gradient2)
    # High(正值)=红色, Low(负值)=蓝色, Mid(0)=白色
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", 
      midpoint = 0, 
      limits = limit_range,
      name = "Change vs CTRL"
    ) +
    
    labs(
      title = "Dynamic Changes Relative to Control",
      subtitle = "Red = Enriched in EGF, Blue = Depleted in EGF (Numbers show count difference)",
      x = "Phospho-site Location", 
      y = "Total Protein Location"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

# --- 执行绘图 ---
cat(">>> 正在生成差异热图 (Delta Heatmap)...\n")
p_delta <- plot_delta_grid(full_bias_df)

# 显示
print(p_delta)

# 保存
ggsave("Fig_Spatial_Bias_Delta_Heatmap.png", p_delta, width = 10, height = 4, dpi = 300)
cat(">>> 图片已保存：Fig_Spatial_Bias_Delta_Heatmap.png\n")

library(tidyverse)
library(ggplot2)

# ==============================================================================
# 终极直观方案：动态迁移气泡图 (Major Shifts Only)
# ==============================================================================
library(tidyverse)
library(ggplot2)
library(ggalluvial)

plot_time_layered_sankey <- function(bias_df) {
  
  # 1. 辅助函数：三大阵营归类
  get_group <- function(fr) {
    case_when(
      fr %in% c("FR1", "FR2") ~ "Cytoplasm",
      fr %in% c("FR3", "FR4") ~ "Membrane",
      fr %in% c("FR5", "FR6") ~ "Nucleus"
    )
  }
  
  # 2. 数据准备：计算增量 (Delta > 0)
  # 统计每个时间点、每个流向的数量
  base_counts <- bias_df %>%
    mutate(
      Group_Prot = get_group(Main_Loc_Prot),
      Group_Phos = get_group(Main_Loc_Phos)
    ) %>%
    group_by(Timepoint, Group_Prot, Group_Phos) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(Timepoint, Group_Prot, Group_Phos, fill = list(Count = 0))
  
  # 提取 Control 作为基准
  ctrl_counts <- base_counts %>%
    filter(Timepoint == "EGF_CTRL") %>%
    select(Group_Prot, Group_Phos, Count_Ctrl = Count)
  
  # 计算 Delta，并准备绘图数据
  sankey_data <- base_counts %>%
    filter(Timepoint != "EGF_CTRL") %>%
    left_join(ctrl_counts, by = c("Group_Prot", "Group_Phos")) %>%
    mutate(Delta = Count - Count_Ctrl) %>%
    # 过滤：只保留发生迁移且数量增加的事件
    filter(Group_Prot != Group_Phos) %>%
    filter(Delta > 0) %>%
    # 确保时间顺序
    mutate(Timepoint = factor(Timepoint, levels = c("EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min")))
  
  # 3. 绘图
  # axis1 = 起点 (总蛋白), axis2 = 终点 (磷酸化)
  # y = 数量 (粗细)
  # fill = 时间 (颜色)
  
  p <- ggplot(sankey_data,
              aes(y = Delta, 
                  axis1 = Group_Prot, 
                  axis2 = Group_Phos)) +
    
    # A. 画流动的彩色丝带
    # 这里的关键是把 Timepoint 映射给 fill
    geom_alluvium(aes(fill = Timepoint), 
                  width = 1/5, 
                  alpha = 0.8,       # 稍微透明一点，防止遮挡
                  curve_type = "sigmoid", # S形曲线，更顺滑
                  color = "white", size = 0.2) + # 丝带间加白色描边，区分度更高
    
    # B. 画两边的柱子 (节点)
    geom_stratum(width = 1/5, fill = "grey90", color = "grey30") +
    
    # C. 加节点文字
    geom_text(stat = "stratum", 
              aes(label = after_stat(stratum)), 
              size = 4, fontface = "bold") +
    
    # D. 颜色设置 (鲜明的对比色)
    scale_fill_manual(values = c(
      "EGF_2min" = "#FFD700",  # 金黄 (早期)
      "EGF_8min" = "#FF7F00",  # 橙色
      "EGF_20min" = "#E41A1C", # 鲜红 (中期爆发)
      "EGF_90min" = "#377EB8"  # 蓝色 (晚期转录)
    ), name = "Stimulation Time") +
    
    # E. 外观调整
    scale_x_discrete(limits = c("Total Protein\n(Origin)", "Phospho-site\n(Destination)"), expand = c(.1, .1)) +
    labs(
      title = "Time-Resolved Spatial Translocation",
      subtitle = "Colored streams represent the number of induced shifts at each timepoint",
      y = "Number of Induced Sites"
    ) +
    theme_void() + # 去掉背景
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
      axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
      legend.position = "bottom",
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

# --- 执行与保存 ---
cat(">>> 正在生成时间分层桑基图...\n")
p_layered <- plot_time_layered_sankey(full_bias_df)

print(p_layered)
ggsave("Fig_Sankey_Time_Layered.png", p_layered, width = 9, height = 7, dpi = 300)
cat(">>> 图片已保存：Fig_Sankey_Time_Layered.png\n")



plot_bubble_final <- function(bias_df) {
  
  # --- 1. 数据准备 (3组分, Delta > 0) ---
  get_group <- function(fr) {
    case_when(
      fr %in% c("FR1", "FR2") ~ "Cytoplasm",
      fr %in% c("FR3", "FR4") ~ "Membrane",
      fr %in% c("FR5", "FR6") ~ "Nucleus"
    )
  }
  
  base_counts <- bias_df %>%
    mutate(
      Group_Prot = get_group(Main_Loc_Prot),
      Group_Phos = get_group(Main_Loc_Phos)
    ) %>%
    group_by(Timepoint, Group_Prot, Group_Phos) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(Timepoint, Group_Prot, Group_Phos, fill = list(Count = 0))
  
  ctrl_counts <- base_counts %>%
    filter(Timepoint == "EGF_CTRL") %>%
    select(Group_Prot, Group_Phos, Count_Ctrl = Count)
  
  plot_data <- base_counts %>%
    filter(Timepoint != "EGF_CTRL") %>%
    left_join(ctrl_counts, by = c("Group_Prot", "Group_Phos")) %>%
    mutate(
      Delta = Count - Count_Ctrl,
      Shift_Type = paste(Group_Prot, "→", Group_Phos) # 这里的箭头是文本标签
    ) %>%
    filter(Group_Prot != Group_Phos) %>%
    filter(Delta > 0) %>%
    mutate(Timepoint = factor(Timepoint, levels = c("EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min")))
  
  # --- 2. 绘图 ---
  p <- ggplot(plot_data, aes(x = Timepoint, y = Shift_Type)) +
    
    # A. 气泡 (放大版)
    geom_point(aes(size = Delta, color = Delta), alpha = 0.9) +
    
    # B. 数字 (居中，加粗)
    geom_text(aes(label = Delta), 
              color = "black", 
              size = 4,         # 字号
              fontface = "bold", 
              vjust = 0.5) +    # 垂直居中
    
    # C. 颜色设置 (暖色调：浅黄 -> 深红)
    scale_color_gradient(low = "#FFD54F", high = "#D32F2F") +
    
    # D. 尺寸设置 (关键：设置最小尺寸为12，保证数字放得下)
    scale_size(range = c(12, 22), guide = "none") + # 去掉Size图例，因为图上有数字
    
    # E. 极简外观
    labs(x = NULL, y = NULL) + # 去掉轴标题，节省空间
    theme_bw() +
    theme(
      # 坐标轴文字加大加粗
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      
      # 去掉背景网格和刻度线
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "grey80", size = 1),
      
      # 图例设置
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

# --- 执行与保存 ---
cat(">>> 正在生成最终版气泡图...\n")
p_bubble_big <- plot_bubble_final(full_bias_df)

print(p_bubble_big)
ggsave("Fig_Translocation_Bubble_Final1.png", p_bubble_big, width = 8, height = 5, dpi = 300)
cat(">>> 图片已保存：Fig_Translocation_Bubble_Final.png\n")

# ==============================================================================
# 2. 磷酸化动态追踪：真迁移 (True Translocation) 与 桑基图
# ==============================================================================
library(ggalluvial)
library(dplyr)
library(ggplot2)

# --- A. 构建全时间点追踪矩阵 (ID Tracking Matrix) ---
# 前提：确保 get_main_location 函数已定义且 df_phos 在内存中
timepoints <- c("EGF_CTRL", "EGF_2min", "EGF_8min", "EGF_20min") 

# 1. 循环计算每个时间点的位置
track_list <- list()

cat(">>> 正在为每个时间点定位 ID...\n")
for (tp in timepoints) {
  loc_res <- get_main_location(df_phos, tp) 
  
  if(!is.null(loc_res)) {
    simp <- loc_res %>% 
      dplyr::select(ID, Main_Loc) %>%
      dplyr::rename(!!paste0("Loc_", tp) := Main_Loc)
    
    track_list[[tp]] <- simp
  }
}

# 2. 合并数据 (逻辑修正：使用 full_join 防止数据因单点缺失被误删)
# 这样即使某位点在 20min 没测到，只要 CTRL 和 8min 在，依然可以分析
full_track_df <- Reduce(function(x, y) full_join(x, y, by = "ID"), track_list)

# 提取 Gene Name
full_track_df$Gene <- sub("_.*", "", full_track_df$ID)

cat(">>> 追踪矩阵构建完成，共追踪位点数:", nrow(full_track_df), "\n")

# --- B. 筛选“真正的迁移位点” (Filtering for Shifters) ---
# 修改这里可以分析不同时间段，例如 CTRL vs 8min
t_start <- "Loc_EGF_CTRL"
t_end   <- "Loc_EGF_8min"

shifters <- full_track_df %>%
  filter(.data[[t_start]] != .data[[t_end]]) %>% # 位置不同
  filter(!is.na(.data[[t_start]]) & !is.na(.data[[t_end]])) # 且两个时间点都有有效值

cat(">>> 在", t_start, "到", t_end, "期间，发生真迁移的位点数:", nrow(shifters), "\n")

# --- C. 桑基图可视化 (Sankey Plot) ---
sankey_data <- shifters %>%
  group_by(.data[[t_start]], .data[[t_end]]) %>% 
  summarise(Count = n(), .groups = "drop")

# 重命名列以便绘图
names(sankey_data)[names(sankey_data) == t_start] <- "Start"
names(sankey_data)[names(sankey_data) == t_end]   <- "End"

# 简化区室名称 (请根据你的实验设计确认 FR 对应关系)
simplify_loc <- function(x) {
  case_when(
    x %in% c("FR1", "FR2") ~ "Cytoplasm",
    x %in% c("FR3", "FR4") ~ "Membrane",
    x %in% c("FR5", "FR6") ~ "Nucleus",
    TRUE ~ "Other"
  )
}

sankey_data <- sankey_data %>%
  mutate(Start_Group = simplify_loc(Start),
         End_Group = simplify_loc(End))

p_sankey <- ggplot(sankey_data,
                   aes(y = Count, axis1 = Start_Group, axis2 = End_Group)) +
  geom_alluvium(aes(fill = Start_Group), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Start (CTRL)", "End (8min)"), expand = c(.05, .05)) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "True Translocation: Site-Specific Tracking",
       subtitle = paste("Tracking Phospho-sites from", t_start, "to", t_end),
       y = "Number of Phospho-sites") +
  theme_void()

print(p_sankey)
ggsave("Fig_Phos_True_Translocation_Corrected.png", p_sankey, width = 8, height = 6)

# --- D. 蛋白归属统计 (Gene Counting) ---
gene_stats <- shifters %>%
  count(Gene, sort = TRUE, name = "Count") 

cat(">>> 发生迁移位点最多的 Top 10 蛋白:\n")
print(head(gene_stats, 10))

p_gene <- ggplot(head(gene_stats, 15), aes(x = reorder(Gene, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top Proteins with Translocating Phospho-sites",
       x = "Gene", y = "Count of Shifting Sites") +
  theme_bw()

print(p_gene)





library(dplyr)
library(tidyverse)

# ==============================================================================
# 修复 Step A: 重新构建全追踪矩阵 (必须显式包含 90min)
# ==============================================================================

# 1. 定义完整的时间点列表 (千万别漏了 CTRL 和 90min)
# 这里的名字必须和你 get_main_location 函数里能识别的名字一致
all_timepoints <- c("EGF_CTRL", "EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min")

track_list <- list()

cat(">>> [修复模式] 正在重新抓取所有时间点的数据 (含 90min)...\n")

for (tp in all_timepoints) {
  # 调用你之前的定位函数
  # 注意：这里假设你的 df_phos 在内存里
  loc_res <- get_main_location(df_phos, tp) 
  
  if(!is.null(loc_res)) {
    # 提取并重命名列
    simp <- loc_res %>% 
      dplyr::select(ID, Main_Loc) %>%
      dplyr::rename(!!paste0("Loc_", tp) := Main_Loc) # 动态命名，如 Loc_EGF_90min
    
    track_list[[tp]] <- simp
    cat(sprintf("   - 成功抓取: %s (行数: %d)\n", tp, nrow(simp)))
  } else {
    cat(sprintf("   - [警告] 无法抓取: %s (函数返回 NULL)\n", tp))
  }
}

# 2. 合并数据 (使用 full_join 防止数据丢失)
cat(">>> 正在合并数据矩阵...\n")
full_track_df <- Reduce(function(x, y) full_join(x, y, by = "ID"), track_list)

# 提取基因名
full_track_df$Gene <- sub("_.*", "", full_track_df$ID)

# 3. 最终核查
if ("Loc_EGF_90min" %in% colnames(full_track_df)) {
  cat("\n>>> [成功] 'Loc_EGF_90min' 列现在存在了！\n")
  cat(sprintf(">>> 全矩阵总行数: %d\n", nrow(full_track_df)))
} else {
  stop(">>> [失败] 依然没有 90min 列。请检查 get_main_location('EGF_90min') 是否能正常工作。")
}

library(ggalluvial)

# ==============================================================================
# 修复 Step B & C: 数据清洗与绘图 (Visual Upgrade)
# ==============================================================================

# --- 1. 准备绘图数据 (包含区室合并逻辑) ---
cat(">>> 正在准备绘图数据...\n")

baseline <- "EGF_CTRL"
# 要展示的受激时间点
stim_times <- c("EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min") 

# 定义区室合并函数
simplify_loc <- function(x) {
  case_when(
    grepl("FR1|FR2", x) ~ "Cytoplasm",
    grepl("FR3|FR4", x) ~ "Membrane",
    grepl("FR5|FR6", x) ~ "Nucleus",
    TRUE ~ "Other" 
  )
}

plot_data_list <- list()

for (tp in stim_times) {
  col_base <- paste0("Loc_", baseline)
  col_stim <- paste0("Loc_", tp)
  
  if (!col_stim %in% colnames(full_track_df)) next # 跳过缺失列
  
  temp_df <- full_track_df %>%
    dplyr::select(ID, Gene, Start_Raw = !!sym(col_base), End_Raw = !!sym(col_stim)) %>%
    filter(!is.na(Start_Raw) & !is.na(End_Raw)) %>%
    filter(Start_Raw != End_Raw) %>% # 原始位置必须不同
    mutate(
      Start = simplify_loc(Start_Raw), # 合并区室
      End   = simplify_loc(End_Raw)
    ) %>%
    filter(Start != End) %>% # 合并后位置也必须不同 (防止 FR1->FR2 这种假阳性)
    mutate(Timepoint = tp)
  
  plot_data_list[[tp]] <- temp_df
}

sankey_data_final <- bind_rows(plot_data_list)
sankey_data_final <- sankey_data_final %>%
  mutate(Group_Prot = Start, Group_Phos = End, Delta = 1)

cat(sprintf(">>> 绘图数据准备就绪。共包含 %d 个跨区室迁移事件。\n", nrow(sankey_data_final)))

# --- 2. 绘图函数 (高对比度 + 智能透明度) ---
plot_time_layered_sankey_v2 <- function(sankey_data) {
  
  # 确保顺序
  sankey_data$Timepoint <- factor(sankey_data$Timepoint, 
                                  levels = c("EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min"))
  
  # 高对比配色
  custom_colors <- c(
    "EGF_2min"  = "#FFC107", # 黄 (亮)
    "EGF_8min"  = "#FF5722", # 橙
    "EGF_20min" = "#C2185B", # 红 (深)
    "EGF_90min" = "#004D40"  # 蓝绿 (极深，高辨识度)
  )
  
  p <- ggplot(sankey_data,
              aes(y = Delta, axis1 = Group_Prot, axis2 = Group_Phos)) +
    geom_alluvium(aes(fill = Timepoint, alpha = Timepoint), 
                  width = 1/6, curve_type = "sigmoid", color = NA) + 
    geom_stratum(width = 1/6, fill = "#F0F0F0", color = "#666666", size=0.3) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
    scale_fill_manual(values = custom_colors, name = "Stimulation Phase") +
    # 智能透明度：90min最不透明
    scale_alpha_manual(values = c("EGF_2min"=0.5, "EGF_8min"=0.6, "EGF_20min"=0.7, "EGF_90min"=0.95), guide="none") +
    scale_x_discrete(limits = c("Start Location\n(CTRL)", "New Location\n(Stimulated)"), expand = c(.08, .08)) +
    labs(title = "Dynamics of Phospho-site Translocation", 
         subtitle = "Tracking cross-compartment shifts (Merged Regions)", y = NULL) +
    theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(hjust=0.5, face="bold"),
          plot.margin = margin(20, 20, 20, 20))
  
  return(p)
}

# --- 3. 出图 ---
p_final <- plot_time_layered_sankey_v2(sankey_data_final)
print(p_final)
ggsave("Fig_Phos_Dynamics_Final_Merged.png", p_final, width = 9, height = 6, dpi = 300, bg="white")
cat(">>> 图片已保存：Fig_Phos_Dynamics_Final_Merged.png\n")