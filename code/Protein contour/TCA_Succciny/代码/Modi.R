# --- [Cell 1] Load & Evaluate Strict Subset ---

library(data.table)
library(dplyr)

# 1. 读取你之前代码保存的文件
# 注意：请确保文件路径正确，如果不在当前目录，需要 setwd()
target_file <- "TaskA_Strict_Pairs_Subset.csv"

if(file.exists(target_file)) {
  strict_data <- read.csv(target_file, stringsAsFactors = FALSE)
  
  cat("✅ 成功读取:", target_file, "\n")
  cat("--- 数据维度 ---\n")
  print(dim(strict_data))
  
  cat("\n--- 列名检查 ---\n")
  print(colnames(strict_data))
  
  # 2. 关键检查：既然之前已经清洗过，这里应该有 'Clean_Seq' 列
  if("Clean_Seq" %in% colnames(strict_data)) {
    cat("\n✅ 发现已清洗的序列列: 'Clean_Seq'\n")
    print(head(strict_data$Clean_Seq))
  } else {
    cat("\n⚠️ 警告: 没找到 'Clean_Seq' 列，可能需要重新清洗。\n")
  }
  
  # 3. 关键检查：各区室的样本分布
  cat("\n--- 各区室样本量分布 (Sample Size per Compartment) ---\n")
  # 这里使用的是修饰后的位置 Predicted_Location
  loc_counts <- table(strict_data$Predicted_Location)
  print(loc_counts)
  
  # 4. 决策逻辑
  if(nrow(strict_data) < 100) {
    cat("\n🛑 导师建议: 样本量少于 100，做 PCA 效果可能极差。建议退回使用 'high_conf_results'。\n")
  } else {
    cat("\n🟢 导师建议: 样本量尚可，可以尝试进行 PCA 分析！\n")
  }
  
} else {
  cat("❌ 错误: 找不到文件", target_file, "。请确认你之前的代码是否成功保存了该文件。")
}

# --- [Cell 2] Strict Subset: Organelle-Specific PCA ---

library(Peptides)
library(dplyr)
library(ggplot2)

# 1. 筛选样本量充足的 Top 5 区室
# -----------------------------------------------------------
# 基于刚才 Cell 1 的输出，我们知道只有这几个够数
target_locs <- c("Ribosome", "Nucleus", "ER", "Mitochondria", "Cytosol")

filtered_strict <- strict_data %>%
  filter(Predicted_Location %in% target_locs)

cat("--- 筛选后用于 PCA 的样本量 ---\n")
print(table(filtered_strict$Predicted_Location))

# 2. 计算理化性质 (Feature Calculation)
# -----------------------------------------------------------
cat("\n正在计算理化性质...\n")

# 电荷 (Charge at pH 7.4)
filtered_strict$Charge <- charge(filtered_strict$Clean_Seq, pH = 7.4)
# 疏水性 (Hydrophobicity - Kyte-Doolittle scale)
filtered_strict$Hydrophobicity <- hydrophobicity(filtered_strict$Clean_Seq, scale = "KyteDoolittle")
# 分子量 (Molecular Weight)
filtered_strict$MW <- mw(filtered_strict$Clean_Seq)
# 脂肪族指数 (Aliphatic Index - 热稳定性指标)
filtered_strict$Aliphatic <- aIndex(filtered_strict$Clean_Seq)
# 序列长度
filtered_strict$Length <- nchar(filtered_strict$Clean_Seq)

# 3. 运行 PCA (Principal Component Analysis)
# -----------------------------------------------------------
# 只选取数值列进行分析
pca_input <- filtered_strict %>% 
  select(Charge, Hydrophobicity, MW, Aliphatic, Length)

# scale. = TRUE 至关重要，因为 MW (几千) 和 Charge (几) 量级不同
pca_res <- prcomp(pca_input, scale. = TRUE)
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

# 4. 可视化 (Visualization)
# -----------------------------------------------------------
pca_df <- as.data.frame(pca_res$x)
pca_df$Location <- filtered_strict$Predicted_Location

p_strict_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Location)) +
  geom_point(size = 3, alpha = 0.7) +
  # 加个圈 (Ellipse) 帮助看清楚聚类，置信区间 0.8
  stat_ellipse(level = 0.8, type = "norm") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Physicochemical Landscape of Translocated Peptides",
       subtitle = "Do modifications in Mitochondria look different from Nucleus?",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"))

print(p_strict_pca)
ggsave("Figure_Strict_Subset_PCA.tiff", p_strict_pca, width = 8, height = 6)

cat("✅ [Cell 2] 完成。请查看生成的 PCA 图。\n")



# --- [Cell 3] Motif Analysis (C-terminal Alignment) --

# 1. 检查并加载包


# 确保 filtered_strict 存在 (来自 Cell 2)
if(!exists("filtered_strict")) {
  stop("❌ 请先运行 Cell 2 生成 filtered_strict 数据")
}

# 2. 序列截取 (Extract C-terminal 7-mer)
# 我们取最后 7 个氨基酸。如果肽段短于 7，就丢掉。
motif_len <- 7

motif_prep <- filtered_strict %>%
  filter(nchar(Clean_Seq) >= motif_len) %>%
  mutate(C_Term_Seq = str_sub(Clean_Seq, -motif_len, -1))

cat("--- C端序列预览 (前5个) ---\n")
print(head(motif_prep$C_Term_Seq))

# 3. 提取对比组：线粒体 vs 细胞核
# -----------------------------------------------------------
seqs_mito <- motif_prep %>% filter(Predicted_Location == "Mitochondria") %>% pull(C_Term_Seq)
seqs_nuc <- motif_prep %>% filter(Predicted_Location == "Nucleus") %>% pull(C_Term_Seq)

cat(paste("线粒体序列数:", length(seqs_mito), "\n"))
cat(paste("细胞核序列数:", length(seqs_nuc), "\n"))

# 4. 绘图 (Sequence Logo)
# -----------------------------------------------------------
# 只有当两组样本都足够多 (>10) 时才画
if(length(seqs_mito) > 10 & length(seqs_nuc) > 10) {
  
  gg_list <- list(
    "Mitochondria" = seqs_mito,
    "Nucleus" = seqs_nuc
  )
  
  p_motif <- ggseqlogo(gg_list, ncol = 1) +
    theme_bw() +
    labs(title = "C-Terminal Motifs Comparison", 
         subtitle = "Do organelles have specific recognition motifs?")
  
  print(p_motif)
  ggsave("Figure_Motif_Mito_vs_Nucleus.tiff", p_motif, width = 6, height = 6)
  cat("✅ [Cell 3] Motif 分析完成！请查看生成的 Logo 图。\n")
  
} else {
  cat("⚠️ [Cell 3] 警告: 某一组序列少于10条，无法生成 Motif 图。\n")
}


# --- [Cell 4] The Final Showdown: Global Mod vs Unmod ---
set.seed(42)

# 1. 准备红队 (Modified Group - Global)
# -----------------------------------------------------------
# 必须使用 high_conf_results (全集)，而不是 strict_data (子集)
if(!exists("high_conf_results")) {
  stop("❌ 错误: 内存中找不到 'high_conf_results'。")
}

# 确保有纯序列列
if(!"Clean_Seq" %in% colnames(high_conf_results)) {
  high_conf_results$Clean_Seq <- str_remove(high_conf_results$Sequence, ":\\s*\\d+$")
}

group_mod <- high_conf_results %>%
  select(Clean_Seq) %>%
  mutate(Type = "Modified")

# 2. 准备蓝队 (Unmod Group - Global Background)
# -----------------------------------------------------------
if(!exists("unmod_data")) {
  if(file.exists("unmod_pep.rds")) {
    unmod_data <- readRDS("unmod_pep.rds")
  } else {
    stop("❌ 找不到 unmod_pep.rds")
  }
}

# 提取并清洗背景序列
all_unmod_seqs <- rownames(fData(unmod_data))
all_unmod_clean <- str_remove(all_unmod_seqs, ":\\s*\\d+$")

# 剔除掉那些其实是修饰的 (Clean Background)
pure_unmod_pool <- setdiff(all_unmod_clean, group_mod$Clean_Seq)

# 随机抽样 (Downsampling)
n_sample <- nrow(group_mod)
# 防止背景组不够抽
if(length(pure_unmod_pool) >= n_sample) {
  sampled_unmod <- sample(pure_unmod_pool, n_sample)
} else {
  sampled_unmod <- pure_unmod_pool
}

group_unmod <- data.frame(
  Clean_Seq = sampled_unmod,
  Type = "Unmodified",
  stringsAsFactors = FALSE
)

# 3. 合并与计算
# -----------------------------------------------------------
df_global <- bind_rows(group_mod, group_unmod)
cat("--- 全局对比样本量 (Global N) ---\n")
print(table(df_global$Type))

cat("\n正在计算全局数据的理化性质 (N =", nrow(df_global), ")... \n")
# 这里也要重新计算一遍，确保无误
df_global$Hydrophobicity <- hydrophobicity(df_global$Clean_Seq, scale = "KyteDoolittle")
df_global$Charge <- charge(df_global$Clean_Seq, pH = 7.4)
df_global$Length <- nchar(df_global$Clean_Seq)

# 4. 关键结果：统计检验与箱线图
# -----------------------------------------------------------
# Wilcoxon 检验
test_hydro <- wilcox.test(Hydrophobicity ~ Type, data = df_global)
cat("\n📢 [结论] 疏水性差异 P-value:", test_hydro$p.value, "\n")

if(test_hydro$p.value < 0.05) {
  cat("✅ 显著！修饰组与未修饰组在疏水性上有显著差异。\n")
} else {
  cat("ℹ️ 不显著。两组疏水性类似。\n")
}

# 绘制箱线图
p_box <- ggplot(df_global, aes(x = Type, y = Hydrophobicity, fill = Type)) +
  geom_boxplot(outlier.alpha = 0.05, outlier.size = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("Modified" = "#E41A1C", "Unmodified" = "#999999")) +
  labs(title = "Global Comparison: Hydrophobicity", 
       subtitle = paste0("Wilcoxon P-value: ", format.pval(test_hydro$p.value)),
       y = "Hydrophobicity Index (Kyte-Doolittle)")

print(p_box)
ggsave("Figure_Global_Boxplot_Hydro.tiff", p_box, width = 6, height = 6)

cat("✅ [Cell 4] 全部完成！请查看箱线图。\n")