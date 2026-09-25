# 设置你的工作目录
setwd("D:\\博士\\TCA_Succciny\\代码")

# --- 2. 加载数据 ---
# 读取琥珀酰化数据 RDS
raw_data <- readRDS("succ_pep.rds")
print(paste("成功加载数据，共", nrow(raw_data), "条肽段。"))

# --- 3. 质量控制 PCA (QC Check) ---
# 这一步是为了告诉师兄：我的数据质量没问题，Marker 分开了。
# 基于 TMT Intensity
qc_pca <- plot2D(raw_data, fcol = "markers", method = "PCA", plot = FALSE)
# 我们自己画好看一点
qc_df <- data.frame(PC1 = qc_pca[,1], PC2 = qc_pca[,2], Marker = fData(raw_data)$markers)
# 只画已知 Marker
qc_df_marker <- qc_df %>% filter(Marker != "unknown")

p_qc <- ggplot(qc_df_marker, aes(x=PC1, y=PC2, color=Marker)) +
  geom_point(alpha=0.7) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "QC: Marker Protein Distribution (Based on Intensity)")
print(p_qc)
ggsave("Step1_QC_PCA.pdf", p_qc, width = 8, height = 6)


# --- 4. 核心任务：SVM 预测定位 ---
# 我们使用 pRoloc 的 SVM 进行预测
cat("正在进行 SVM 预测...\n")
# 注意：如果 cost/sigma 需要调整，请修改这里
svm_model <- svmClassification(raw_data, cost = 16, sigma = 0.1)

# 提取预测结果
# 这里的 threshold 设为 0.7，保证我们分析的都是“靠谱”的
predictions <- fData(svm_model)$svm
scores <- fData(svm_model)$svm.scores

# 构建结果大表
result_df <- data.frame(
  Sequence_Raw = rownames(fData(svm_model)),
  Predicted_Location = predictions,
  Score = scores,
  stringsAsFactors = FALSE
)

# 筛选高置信度结果 (这是我们后续分析的基础！)
high_conf_df <- result_df %>% filter(Score > 0.7)

cat(paste("SVM 预测完成。原始数据:", nrow(result_df), 
          "| 高置信度数据:", nrow(high_conf_df), "\n"))
print(table(high_conf_df$Predicted_Location))


# --- 5. 师兄的任务：理化性质分析 (Physicochemical Analysis) ---

cat("\n开始特征提取...\n")

# 5.1 序列清洗：去掉修饰标记 (如 K123) 变成纯字母序列
high_conf_df$Clean_Seq <- str_remove(high_conf_df$Sequence_Raw, "[:\\.].*$") 

# 5.2 计算 5 大特征
# 使用 Peptides 包
high_conf_df$Charge <- charge(high_conf_df$Clean_Seq, pH = 7.4)
high_conf_df$Hydrophobicity <- hydrophobicity(high_conf_df$Clean_Seq, scale = "KyteDoolittle")
high_conf_df$MW <- mw(high_conf_df$Clean_Seq)
high_conf_df$Aliphatic <- aIndex(high_conf_df$Clean_Seq)
high_conf_df$pI <- pI(high_conf_df$Clean_Seq)
high_conf_df$Length <- nchar(high_conf_df$Clean_Seq)
high_conf_df$Type <- "Modified" # 打上标签

# 只选取样本量充足的细胞器 (>50)
target_locs <- names(which(table(high_conf_df$Predicted_Location) > 50))
pca_input_df <- high_conf_df %>% filter(Predicted_Location %in% target_locs)

# --- 5.3 宏观视角：PCA 分析 (Multivariate) ---
cat("运行 PCA...\n")
pca_matrix <- pca_input_df %>% select(Charge, Hydrophobicity, MW, Aliphatic, pI, Length)
pca_res_chem <- prcomp(pca_matrix, scale. = TRUE)

pca_plot_df <- data.frame(
  PC1 = pca_res_chem$x[,1], 
  PC2 = pca_res_chem$x[,2], 
  Location = pca_input_df$Predicted_Location
)

var_exp <- round(summary(pca_res_chem)$importance[2, 1:2] * 100, 1)

p_chem_pca <- ggplot(pca_plot_df, aes(x=PC1, y=PC2, color=Location)) +
  geom_point(alpha=0.6, size=2) +
  stat_ellipse(level = 0.8) + 
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Physicochemical Landscape (Global PCA)",
       subtitle = "If overlapping, individual features might still differ.",
       x = paste0("PC1 - ", var_exp[1], "%"),
       y = paste0("PC2 - ", var_exp[2], "%"))

print(p_chem_pca)
ggsave("Step2_Physicochem_PCA.pdf", p_chem_pca, width = 8, height = 6)


# --- 5.4 微观视角：单变量箱线图 (Univariate Boxplots) ---
# 这里回答你的疑问：一个一个看！
cat("运行单变量箱线图分析...\n")

features_to_plot <- c("Charge", "Hydrophobicity", "Aliphatic", "pI")

for (feat in features_to_plot) {
  # 动态画图 - 修复 aes_string 警告，使用 .data[[feat]]
  p_box <- ggplot(pca_input_df, aes(x = Predicted_Location, y = .data[[feat]], fill = Predicted_Location)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    # 添加统计显著性 (Kruskal-Wallis 检验)
    stat_compare_means(label.y.npc = "top", label.x.npc = "left") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    scale_fill_brewer(palette = "Set3") +
    labs(title = paste0("Feature Distribution: ", feat),
         subtitle = "Comparison across organelles",
         y = feat, x = "")
  
  print(p_box)
  ggsave(paste0("Step2_Boxplot_", feat, ".pdf"), p_box, width = 6, height = 5)
}
cat("单变量分析完成。请检查 Step2_Boxplot 系列图片，寻找 P值显著的结果。\n")


# ==============================================================================
# [新增模块] 5B. Modified vs Unmodified 对比 (Global Analysis)
# ==============================================================================
cat("\n正在准备 Unmodified 数据进行对比...\n")

if(file.exists("unmod_pep.rds")) {
  unmod_data <- readRDS("unmod_pep.rds")
  
  # 1. 提取未修饰序列
  unmod_seqs_raw <- rownames(fData(unmod_data))
  unmod_seqs_clean <- str_remove(unmod_seqs_raw, "[:\\.].*$")
  
  # 2. 剔除那些其实是修饰过的
  unmod_seqs_clean <- setdiff(unmod_seqs_clean, high_conf_df$Clean_Seq)
  
  # 3. 随机抽样
  set.seed(123)
  n_sample <- nrow(high_conf_df)
  
  if(length(unmod_seqs_clean) >= n_sample) {
    unmod_sample_seqs <- sample(unmod_seqs_clean, n_sample)
  } else {
    unmod_sample_seqs <- unmod_seqs_clean
  }
  
  # 4. 构建 Unmod 数据框
  unmod_df <- data.frame(
    Clean_Seq = unmod_sample_seqs,
    Type = "Unmodified",
    stringsAsFactors = FALSE
  )
  
  # 5. 计算 Unmod 特征
  unmod_df$Charge <- charge(unmod_df$Clean_Seq, pH = 7.4)
  unmod_df$Hydrophobicity <- hydrophobicity(unmod_df$Clean_Seq, scale = "KyteDoolittle")
  unmod_df$MW <- mw(unmod_df$Clean_Seq)
  unmod_df$Aliphatic <- aIndex(unmod_df$Clean_Seq)
  unmod_df$pI <- pI(unmod_df$Clean_Seq)
  unmod_df$Length <- nchar(unmod_df$Clean_Seq)
  
  # 6. 合并红蓝两队
  cols_to_keep <- c("Type", "Charge", "Hydrophobicity", "MW", "Aliphatic", "pI", "Length")
  
  # --- FIX START: 增加安全检查 ---
  # 防止用户没运行上面的 feature 计算步骤
  if(!"Type" %in% colnames(high_conf_df)) high_conf_df$Type <- "Modified"
  
  required_cols <- setdiff(cols_to_keep, "Type")
  missing_cols <- setdiff(required_cols, colnames(high_conf_df))
  if(length(missing_cols) > 0) {
    stop(paste("错误：high_conf_df 缺少以下特征列，请先运行上方的 '5. 特征提取' 代码块:", 
               paste(missing_cols, collapse=", ")))
  }
  # --- FIX END ---
  
  # 使用 bind_rows + select 替代 rbind，更安全
  global_pca_input <- bind_rows(
    high_conf_df %>% select(all_of(cols_to_keep)),
    unmod_df %>% select(all_of(cols_to_keep))
  )
  
  # 7. 运行 Global PCA
  cat("运行 Modified vs Unmodified PCA...\n")
  pca_matrix_global <- global_pca_input %>% select(Charge, Hydrophobicity, MW, Aliphatic, pI, Length)
  pca_matrix_global <- na.omit(pca_matrix_global)
  
  pca_res_global <- prcomp(pca_matrix_global, scale. = TRUE)
  
  pca_plot_global <- data.frame(
    PC1 = pca_res_global$x[,1], 
    PC2 = pca_res_global$x[,2], 
    Type = global_pca_input$Type[1:nrow(pca_matrix_global)]
  )
  
  var_exp_g <- round(summary(pca_res_global)$importance[2, 1:2] * 100, 1)
  
  p_global_pca <- ggplot(pca_plot_global, aes(x=PC1, y=PC2, color=Type)) +
    geom_point(alpha=0.5, size=2) +
    stat_ellipse(level = 0.95) + 
    theme_bw() +
    scale_color_manual(values = c("Modified" = "red", "Unmodified" = "grey")) +
    labs(title = "PCA: Modified vs Unmodified Peptides",
         subtitle = "Distinct physicochemical properties?",
         x = paste0("PC1 - ", var_exp_g[1], "%"),
         y = paste0("PC2 - ", var_exp_g[2], "%"))
  
  print(p_global_pca)
  ggsave("Step2B_Global_Mod_vs_Unmod_PCA.pdf", p_global_pca, width = 8, height = 6)
  
  # 8. 对应的单变量箱线图
  cat("运行 Modified vs Unmodified 箱线图...\n")
  for (feat in c("Charge", "Hydrophobicity", "Aliphatic", "pI")) {
    # 修复 aes_string 警告，使用 .data[[feat]]
    p_box_g <- ggplot(global_pca_input, aes(x = Type, y = .data[[feat]], fill = Type)) +
      geom_boxplot(outlier.size = 0.5) +
      stat_compare_means(label.y.npc = "top") + 
      theme_bw() +
      scale_fill_manual(values = c("Modified" = "red", "Unmodified" = "grey")) +
      labs(title = paste0("Global Comparison: ", feat))
    
    print(p_box_g)
    ggsave(paste0("Step2B_Global_Boxplot_", feat, ".pdf"), p_box_g, width = 5, height = 5)
  }
  
} else {
  cat("警告：找不到 unmod_pep.rds，无法进行 Mod vs Unmod 对比。\n")
}

# ==============================================================================
# 6. 桑基图 (Sankey Plot) - 极简回归版 (美化修复 + 过滤)
# ==============================================================================
# 仅当 unmod 文件存在时尝试
if(file.exists("unmod_pep.rds")) {
  cat("\n正在处理未修饰数据以绘制桑基图...\n")
  
  # 1. 准备 Unmod 数据 (Home Location)
  if(!exists("unmod_model")) {
    unmod_data <- readRDS("unmod_pep.rds")
    unmod_model <- svmClassification(unmod_data, cost = 16, sigma = 0.125)
  }
  
  # 提取 Unmod 的 Accession 和预测结果
  unmod_raw_ids <- fData(unmod_model)$Master.Protein.Accessions
  # 简单清洗：只取分号前的第一个ID
  unmod_clean_ids <- sapply(strsplit(as.character(unmod_raw_ids), ";"), `[`, 1)
  
  unmod_preds <- data.frame(
    Accession = unmod_clean_ids,
    Location = fData(unmod_model)$svm,
    Score = fData(unmod_model)$svm.scores
  ) %>% filter(Score > 0.5) 
  
  # 确定老家 (Home)
  protein_home <- unmod_preds %>%
    group_by(Accession) %>%
    summarise(Home_Location = names(which.max(table(Location))))
  
  # 2. 准备 Mod 数据 (New Location)
  mod_raw_ids <- fData(svm_model)[high_conf_df$Sequence_Raw, "Master.Protein.Accessions"]
  mod_clean_ids <- sapply(strsplit(as.character(mod_raw_ids), ";"), `[`, 1)
  
  sankey_input_prep <- high_conf_df %>%
    mutate(Accession = mod_clean_ids)
  
  # 3. 合并数据
  sankey_input <- sankey_input_prep %>%
    inner_join(protein_home, by = "Accession")
  
  cat(paste("成功合并数据，用于绘图的路径数量:", nrow(sankey_input), "\n"))
  
  # 4. 绘图 (修复颜色和标签 + 过滤)
  sankey_data <- sankey_input %>%
    group_by(Home_Location, Predicted_Location) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    filter(Count > 20) # [关键修复] 过滤掉数量少于 20 的路径，避免图片杂乱
  
  if(nrow(sankey_data) > 0) {
    # 动态扩展调色板 (解决颜色不够用的报错)
    library(RColorBrewer)
    n_colors <- length(unique(sankey_data$Home_Location))
    # 如果类别超过 Set1 的上限 (9)，则使用 colorRampPalette 自动生成更多颜色
    my_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_colors)
    
    p_sankey <- ggplot(sankey_data,
                       aes(y = Count, axis1 = Home_Location, axis2 = Predicted_Location)) +
      geom_alluvium(aes(fill = Home_Location), width = 1/12, alpha = 0.7) + # 增加透明度，看清流向
      geom_stratum(width = 1/12, fill = "grey90", color = "grey30") +
      
      # [智能标签]: 只有当该层级的数量 (Count) 大于总数的 1% 时才显示标签，避免重叠
      # 注意：ggalluvial 的 after_stat(n) 计算的是总数，这里我们简化处理，
      # 直接用 geom_text 配合 stat="stratum"，如果不显示可以用 size 控制
      geom_label(stat = "stratum", aes(label = after_stat(stratum)), 
                 size = 3,  # 调小字号
                 alpha = 0.9) + # 标签背景略微透明
      
      scale_x_discrete(limits = c("Unmodified (Home)", "Succinylated (New)"), expand = c(.05, .05)) +
      scale_fill_manual(values = my_colors) + # 使用自定义扩展色板
      ggtitle(paste0("Protein Translocation Flow (N=", sum(sankey_data$Count), ")")) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.y = element_blank(), # 隐藏 Y 轴数字，让图更干净
            axis.ticks = element_blank(),
            panel.grid = element_blank())
    
    print(p_sankey)
    ggsave("Figure_Translocation_Sankey.TIFF", p_sankey, width = 10, height = 8, dpi = 300)
    cat("✅ 桑基图已保存为 Figure_Translocation_Sankey.TIFF (已修复颜色和重叠问题，并应用了 Count > 20 的过滤)\n")
  } else {
    cat("❌ 桑基图数据仍为空。请检查 Mod 和 Unmod 数据的 Accession 是否真的有重叠。\n")
  }
} else {
  cat("警告：找不到 unmod_pep.rds，无法进行 Mod vs Unmod 对比，跳过桑基图绘制。\n")
}

cat("\n所有分析结束。请检查工作目录下的生成图片。\n")