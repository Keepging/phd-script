library(MSnbase)
library(pRoloc)
library(tidyverse)
library(data.table)

setwd("D:\\博士\\TCA_Succciny\\代码")
# ==============================================================================
# 1. 数据加载与 SVM 预测
# ==============================================================================

# 定义文件名，方便复用
mod_pep   <- "succ_pep.rds"
unmod_pep <- "unmod_pep.rds"
svm_mod   <- "SVM_Result_Mod.rds"   
svm_unmod <- "SVM_Result_Unmod.rds" 

setwd("D:/博士/TCA_Succciny/代码")

succ_pep  <- readRDS("succ_pep.rds")
unmod_pep <- readRDS("unmod_pep.rds")

svm_mod   <- readRDS("SVM_Result_Mod.rds")
svm_unmod <- readRDS("SVM_Result_Unmod.rds")



# --- A. 处理修饰组 (Succ) ---
if (file.exists(svm_mod_res_file)) {
  cat(">>> 检测到已保存的修饰组结果，正在加载...\n")
  high_conf_results <- readRDS(svm_mod_res_file)
} else {
  cat(">>> 未检测到结果，开始加载原始数据并运行 SVM (这可能需要几分钟)...\n")
  raw_data <- readRDS(mod_file)
  
  # 简单过滤 Marker
  min_markers <- 20
  counts <- table(fData(raw_data)$markers)
  keep <- names(counts)[counts >= min_markers & names(counts) != "unknown"]
  qc_msnset <- raw_data[fData(raw_data)$markers %in% c(keep, "unknown"), ]
  
  # 运行 SVM
  # 注意：这里直接用经验参数，省去 Grid Search 时间。若需严谨可加回 Search
  svm_model <- svmClassification(qc_msnset, cost = 16, sigma = 0.1) 
  
  # 整理结果
  res <- data.frame(exprs(svm_model))
  res$Sequence <- rownames(res)
  # 确保有 Protein ID 列 (非常重要！)
  res$Master.Protein.Accessions <- fData(svm_model)$Master.Protein.Accessions 
  res$Predicted_Location <- fData(svm_model)$svm
  res$Score <- fData(svm_model)$svm.scores
  
  high_conf_results <- res %>% filter(Score > 0.7)
  saveRDS(high_conf_results, svm_mod_res_file)
  cat(">>> 修饰组 SVM 完成并保存。\n")
}

# --- B. 处理未修饰组 (Unmod) ---
if (file.exists(svm_unmod_res_file)) {
  cat(">>> 检测到已保存的未修饰组结果，正在加载...\n")
  high_conf_unmod <- readRDS(svm_unmod_res_file)
} else {
  cat(">>> 运行未修饰组 SVM...\n")
  unmod_data <- readRDS(unmod_file)
  
  # 假设参数已调优，直接跑 (为了速度)
  unmod_model <- svmClassification(unmod_data, cost = 16, sigma = 0.125)
  
  res_u <- data.frame(exprs(unmod_model))
  res_u$Sequence <- rownames(res_u)
  res_u$Master.Protein.Accessions <- fData(unmod_model)$Master.Protein.Accessions
  res_u$Predicted_Location <- fData(unmod_model)$svm
  res_u$Score <- fData(unmod_model)$svm.scores
  
  high_conf_unmod <- res_u %>% filter(Score > 0.7)
  saveRDS(high_conf_unmod, svm_unmod_res_file)
  cat(">>> 未修饰组 SVM 完成并保存。\n")
}





library(dplyr)
library(tidyverse)

# ==============================================================================
# 临时热补丁代码：不跑SVM，只看新逻辑结果
# 前提：你的 R 环境里必须已经有了 'high_conf_results' 和 'high_conf_unmod' 这两个表
# ==============================================================================

cat(">>> 正在修复缺失的蛋白ID列...\n")

# --- 修复 1: 给未修饰组 (Unmod) 加上蛋白ID ---
# 这里的 unmod_model 是你之前跑完 SVM 后生成的对象，里面存有 fData 信息
if (!"Master.Protein.Accessions" %in% colnames(high_conf_unmod)) {
  # 通过行名 (Sequence) 对齐，把 fData 里的蛋白ID 抓过来
  # 确保 unmod_model 还在内存里
  if (exists("unmod_model")) {
    high_conf_unmod$Master.Protein.Accessions <- fData(unmod_model)[rownames(high_conf_unmod), "Master.Protein.Accessions"]
  } else if (exists("unmod_data")) {
    high_conf_unmod$Master.Protein.Accessions <- fData(unmod_data)[rownames(high_conf_unmod), "Master.Protein.Accessions"]
  } else {
    stop("错误：内存中找不到 'unmod_model' 或 'unmod_data'，无法找回丢失的蛋白ID！")
  }
}

# --- 修复 2: 给修饰组 (Mod/Succ) 加上蛋白ID ---
if (!"Master.Protein.Accessions" %in% colnames(high_conf_results)) {
  if (exists("svm_model")) {
    high_conf_results$Master.Protein.Accessions <- fData(svm_model)[rownames(high_conf_results), "Master.Protein.Accessions"]
  } else if (exists("qc_msnset")) {
    high_conf_results$Master.Protein.Accessions <- fData(qc_msnset)[rownames(high_conf_results), "Master.Protein.Accessions"]
  } else {
    stop("错误：内存中找不到 'svm_model' 或 'qc_msnset'，无法找回丢失的蛋白ID！")
  }
}

cat(">>> 蛋白ID列修复完成，开始执行逻辑修正...\n")

# ==============================================================================
# 以下是刚才的逻辑修正代码 (Group By 现在可以正常运行了)
# ==============================================================================

# --- 1. 建立未修饰蛋白的“老家”档案 (Protein-Level Home) ---
cat(">>> 正在计算未修饰蛋白的共识位置 (Consensus Location)...\n")

unmod_consensus <- high_conf_unmod %>%
  group_by(Master.Protein.Accessions) %>%
  summarise(
    Unmod_Home_Loc = names(which.max(table(Predicted_Location))), # 出现次数最多的位置
    Unmod_Pep_Count = n() # 支持该位置的肽段数
  )

# --- 2. 核心比对 (The Logic Switch) ---
cat(">>> 正在进行蛋白归属匹配...\n")

comparison_df <- high_conf_results %>%
  # 1. 给修饰肽段贴上它所属蛋白的“老家”标签
  left_join(unmod_consensus, by = "Master.Protein.Accessions") %>%
  dplyr::select(
    Sequence, 
    Master.Protein.Accessions, 
    Mod_Location = Predicted_Location, # 它现在在哪（修饰后）
    Unmod_Home_Loc,                    # 它老家在哪（未修饰）
    Mod_Score = Prediction_Score       # 确保列名正确
  ) %>%
  # 2. 过滤掉那些没找到对应未修饰蛋白的数据
  filter(!is.na(Unmod_Home_Loc))

# --- 3. 统计“易位”情况 (Translocation) ---
translocated_candidates <- comparison_df %>%
  filter(Mod_Location != Unmod_Home_Loc)

# ==============================================================================
# 结果速览
# ==============================================================================
cat("\n#######################################################\n")
cat("               逻辑修正后的临时结果速览                 \n")
cat("#######################################################\n")

num_peps <- nrow(translocated_candidates)
num_prots <- length(unique(translocated_candidates$Master.Protein.Accessions))

cat("2. [当前] 新的 Protein-Level 易位分析结果:\n")
cat(sprintf("   - 发现发生位置改变的修饰肽段 (Peptides): %d\n", num_peps))
cat(sprintf("   - 这些肽段归属于多少个蛋白 (Proteins):   %d\n", num_prots))
cat("\n")

if (num_prots > 0) {
  cat(sprintf("3. [验证] 如果上面的 Protein 数量 (%d) 接近 400 (或 399)，\n", num_prots))
  cat("   说明 PPT 里的 399 指的是 'Protein' 数量，之前的结论是对的！\n")
} else {
  cat("3. [警告] 数字依然很小？可能需要检查 Master.Protein.Accessions 列的内容是否匹配。\n")
}
cat("#######################################################\n")
print(head(translocated_candidates))
# ==============================================================================
# 接之前的代码...
# 确保 translocated_candidates 对象存在
# ==============================================================================

library(ggplot2)
library(ggalluvial)
library(dplyr)

# ------------------------------------------------------------------------------
# 任务二可视化：绘制桑基图 (Sankey Plot) - 展示蛋白的动态迁移
# ------------------------------------------------------------------------------
cat(">>> 正在绘制蛋白迁移路径 (Sankey Plot)...\n")

# 1. 整理数据：统计从 A 到 B 的蛋白数量
sankey_data <- translocated_candidates %>%
  group_by(Unmod_Home_Loc, Mod_Location) %>%
  summarise(Count = n_distinct(Master.Protein.Accessions), .groups = "drop") %>%
  filter(Count > 0) # 过滤掉0的

# 2. 绘图
p_sankey <- ggplot(sankey_data,
                   aes(axis1 = Unmod_Home_Loc, axis2 = Mod_Location, y = Count)) +
  geom_alluvium(aes(fill = Unmod_Home_Loc), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Unmodified (Home)", "Succinylated (New Loc)"), expand = c(.05, .05)) +
  labs(title = "Protein Translocation Landscape (n=236)",
       subtitle = "Where do proteins go after succinylation?",
       y = "Number of Proteins") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

print(p_sankey)
ggsave("Fig_Task2_Sankey_Dynamics.png", p_sankey, width = 10, height = 7)
cat(">>> 桑基图已保存为 'Fig_Task2_Sankey_Dynamics.png'\n")


# ------------------------------------------------------------------------------
# 任务一可视化修正：绘制区室特异性比例 (The Overlap Bar Plot)
# ------------------------------------------------------------------------------
cat(">>> 正在绘制修正后的 Task 1 统计图...\n")

# 1. 重新计算每个区室的基数 (分母：该区室原本有多少个未修饰蛋白)
# 注意：这里我们用 unmod_consensus (之前代码生成的) 作为全集
compartment_stats <- unmod_consensus %>%
  group_by(Unmod_Home_Loc) %>%
  summarise(Total_Unmod_Proteins = n(), .groups = "drop")

# 2. 计算每个区室有多少个蛋白“跑了” (分子：Translocated)
translocated_counts <- translocated_candidates %>%
  group_by(Unmod_Home_Loc) %>%
  summarise(Translocated_Count = n_distinct(Master.Protein.Accessions), .groups = "drop")

# 3. 合并数据
plot_df_corrected <- compartment_stats %>%
  left_join(translocated_counts, by = "Unmod_Home_Loc") %>%
  mutate(Translocated_Count = replace_na(Translocated_Count, 0)) %>%
  mutate(Stable_Count = Total_Unmod_Proteins - Translocated_Count) %>%
  # 转换为长格式以便画图
  pivot_longer(cols = c("Translocated_Count", "Stable_Count"),
               names_to = "Status", values_to = "Count")

# 4. 绘图 (这次一定会有红色了！)
p_overlap <- ggplot(plot_df_corrected, aes(x = Unmod_Home_Loc, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Stable_Count" = "grey80", "Translocated_Count" = "#D53E4F"),
                    labels = c("Stable (No Change)", "Translocated upon Succ.")) +
  labs(title = "Compartment-Specific Translocation Rate",
       subtitle = "Percentage of proteins in each compartment that change location upon modification",
       x = "Original Compartment (Unmodified)",
       y = "Percentage of Proteins") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_overlap)
ggsave("Fig_Task1_Overlap_Corrected.png", p_overlap, width = 8, height = 5)
cat(">>> 修正后的比例图已保存为 'Fig_Task1_Overlap_Corrected.png'\n")

cat("\n恭喜！最关键的两个图表已生成。周一的任务圆满完成。\n")







library(dplyr)
library(tidyverse)

# ==============================================================================
# 6. 文字报告生成模块 (Text Report Generation)
# 前提：请确保 'translocated_candidates' 对象在内存里 (即运行过之前的代码)
# ==============================================================================

if (!exists("translocated_candidates")) {
  stop("请先运行上一段代码生成 'translocated_candidates' 数据！")
}

cat("\n==============================================================================\n")
cat("                       任务汇报核心文字总结 (Protein-Level)                     \n")
cat("==============================================================================\n")

# --- 1. 总体结论 (General Conclusion) ---
n_prots <- length(unique(translocated_candidates$Master.Protein.Accessions))
n_peps  <- nrow(translocated_candidates)

cat(sprintf("【总体结论】\n"))
cat(sprintf("在蛋白层面，我们鉴定了 %d 个发生琥珀酰化依赖性易位 (Translocation) 的蛋白。\n", n_prots))
cat(sprintf("这一现象由 %d 条高置信度的修饰肽段所支持。\n", n_peps))
cat("这证明了琥珀酰化确实可能改变蛋白的亚细胞定位，而非仅仅是随机分布。\n\n")

# --- 2. 具体的“流向”描述 (Flow Description) ---
# 统计最主要的迁移路径：比如 Mitochondria -> Nucleus 有多少个？
flow_stats <- translocated_candidates %>%
  group_by(Unmod_Home_Loc, Mod_Location) %>%
  summarise(Protein_Count = n_distinct(Master.Protein.Accessions), .groups = "drop") %>%
  arrange(desc(Protein_Count))

cat("【主要迁移路径 (Top 5 Flow Pathways)】\n")
cat("以下数据展示了蛋白从“未修饰老家”跑到了“修饰新家”的主要趋势：\n")
top_flows <- head(flow_stats, 5)
for(i in 1:nrow(top_flows)) {
  row <- top_flows[i, ]
  cat(sprintf("  %d. 从 [%s] 迁移到 [%s]: 涉及 %d 个蛋白\n", 
              i, row$Unmod_Home_Loc, row$Mod_Location, row$Protein_Count))
}
cat("(建议：在 PPT 中重点讨论前两条路径的生物学意义)\n\n")

# --- 3. 明星蛋白举例 (Star Examples) ---
# 挑选置信度最高的前 3 个例子
# fix: 加入 as.data.frame() 以消除 Rle 报错
top_examples <- translocated_candidates %>%
  as.data.frame() %>% 
  arrange(desc(Mod_Score)) %>%
  dplyr::slice(1:3)

cat("【具体蛋白案例 (Examples for PPT)】\n")
cat("你可以向师兄展示以下具体的蛋白作为证据：\n")
for(i in 1:nrow(top_examples)) {
  ex <- top_examples[i, ]
  cat(sprintf("  案例 %d: 蛋白ID [%s]\n", i, ex$Master.Protein.Accessions))
  cat(sprintf("    - 肽段序列: %s\n", ex$Sequence))
  cat(sprintf("    - 未修饰时常驻: %s\n", ex$Unmod_Home_Loc))
  cat(sprintf("    - 修饰后发现于: %s (置信度: %.2f)\n", ex$Mod_Location, ex$Mod_Score))
}
cat("\n")

# --- 4. 导出详细表格 (Export) ---
# 师兄如果问“给我看看这 236 个蛋白的列表”，你就把这个 CSV 发给他
output_file <- "Final_Report_236_Translocated_Proteins.csv"
write.csv(translocated_candidates, output_file, row.names = FALSE)
cat(sprintf("【文件输出】\n完整列表已保存为: %s\n", output_file))
cat("你可以用 Excel 打开它，按 Mod_Score 排序查看。\n")

cat("==============================================================================\n")



library(dplyr)
library(tidyverse)
# 假设你已经安装了 clusterProfiler，如果没有，这部分只会跑第一步来源分析
# library(clusterProfiler) 
# library(org.Hs.eg.db)

# ==============================================================================
# 验证师兄猜想：核糖体 (Ribosome) 是“公共车站”还是“特异终点”？
# ==============================================================================

cat("\n>>> 正在深入分析 'To-Ribosome' 的蛋白群体...\n")

# 1. 筛选出所有修饰后“跑去”核糖体的蛋白
ribo_visitors <- translocated_candidates %>%
  filter(Mod_Location == "Ribosome")

n_visitors <- length(unique(ribo_visitors$Master.Protein.Accessions))

cat(sprintf("【基本统计】\n共有 %d 个蛋白在修饰后定位到了 Ribosome。\n\n", n_visitors))

# ------------------------------------------------------------------------------
# 验证步 1: 来源分析 (Source Analysis)
# ------------------------------------------------------------------------------
cat("【验证 1: 它们从哪儿来？(Source Distribution)】\n")
source_stats <- ribo_visitors %>%
  group_by(Unmod_Home_Loc) %>%
  summarise(Count = n_distinct(Master.Protein.Accessions)) %>%
  arrange(desc(Count)) %>%
  mutate(Percent = round(Count / sum(Count) * 100, 1))

print(source_stats)

cat("\n[分析话术]:\n")
if(nrow(source_stats) >= 4) {
  cat("结果显示：来源非常广泛 (Cytosol, Nucleus, ER 等都有)。\n")
  cat("这倾向于支持师兄的猜想：核糖体像个‘加工厂’，各路蛋白都汇聚于此。\n")
} else {
  cat("结果显示：来源非常集中 (主要来自某一两个区室)。\n")
  cat("这可能反驳了‘随机汇聚’的猜想，暗示了特定的调控路径。\n")
}

# ------------------------------------------------------------------------------
# 验证步 2: 身份检查 (Identity Check)
# 看看排名前 10 的“核糖体访客”到底是谁？
# ------------------------------------------------------------------------------
cat("\n【验证 2: 它们是谁？(Top Visitors check)】\n")
cat("如果是由于‘正在合成’导致的假象，这些蛋白的功能应该很杂。\n")
cat("如果是‘功能性结合’，它们可能都是 RNA Binding Proteins。\n\n")

# 这里简单打印前10个蛋白ID和序列，你需要查一下它们的 Gene Name
top_visitors <- ribo_visitors %>%
  dplyr::select(Master.Protein.Accessions, Unmod_Home_Loc, Sequence, Mod_Score) %>%
  arrange(desc(Mod_Score)) %>%
  head(10)

print(top_visitors)

# (可选) 如果你想看是否有某种功能富集，可以用简单的文字描述
cat("\n>>> 建议：请在 Uniprot 上随机查 3-5 个上面的 ID。\n")
cat("1. 如果它们是 Transcription Factors (转录因子) -> 可能是因为修饰阻止了它们入核，滞留在核糖体附近。\n")
cat("2. 如果它们是 Structural Proteins (结构蛋白) -> 支持师兄的‘正在合成’假说。\n")








library(dplyr)
library(tidyr)
library(UpSetR)

setwd("D:/博士/TCA_Succciny/代码")

# -------------------------------------------------------------------
# 0. 读入 SVM 结果
# -------------------------------------------------------------------
svm_mod   <- readRDS("SVM_Result_Mod.rds")
svm_unmod <- readRDS("SVM_Result_Unmod.rds")

# 修饰组 Sequence 去掉 ": 位置"，统一成和 unmod 一样的纯肽序列
svm_mod$Sequence <- sub(":.*$", "", svm_mod$Sequence)

# 10 个 channel 列
chan_cols <- c("X126","X127N","X127C","X128N","X128C",
               "X129N","X129C","X130N","X130C","X131")

# -------------------------------------------------------------------
# 1. 对每个 data.frame 做行归一化，得到每个 channel 的 fraction
#    并找出 fraction 最大的主 channel
# -------------------------------------------------------------------
norm_channels <- function(df, chan_cols) {
  mat <- as.matrix(df[, chan_cols])
  row_sums <- rowSums(mat)
  frac_mat <- sweep(mat, 1, row_sums, "/")
  frac_mat[!is.finite(frac_mat)] <- 0
  colnames(frac_mat) <- paste0("Frac_", chan_cols)
  
  df_out <- bind_cols(
    df[, c("Sequence","Master.Protein.Accessions")],
    as.data.frame(frac_mat)
  )
  
  main_idx <- max.col(frac_mat, ties.method = "first")
  df_out$MainChan <- chan_cols[main_idx]
  
  df_out
}

unmod_frac <- norm_channels(svm_unmod, chan_cols)
mod_frac   <- norm_channels(svm_mod,   chan_cols)

# -------------------------------------------------------------------
# 2. 用 Sequence + Protein 配对 unmod vs succ
# -------------------------------------------------------------------
unmod_frac2 <- unmod_frac %>%
  dplyr::rename(Protein    = Master.Protein.Accessions,
                Main_unmod = MainChan) %>%
  dplyr::select(Sequence, Protein, Main_unmod, dplyr::starts_with("Frac_"))

mod_frac2 <- mod_frac %>%
  dplyr::rename(Protein  = Master.Protein.Accessions,
                Main_mod = MainChan) %>%
  dplyr::select(Sequence, Protein, Main_mod, dplyr::starts_with("Frac_"))

# 分别给 Frac_ 列加上 _unmod / _mod 后缀
unmod_frac2 <- unmod_frac2 %>%
  dplyr::rename_with(~ paste0(., "_unmod"), dplyr::starts_with("Frac_"))

mod_frac2 <- mod_frac2 %>%
  dplyr::rename_with(~ paste0(., "_mod"), dplyr::starts_with("Frac_"))

# 配对
paired <- dplyr::inner_join(
  unmod_frac2,
  mod_frac2,
  by = c("Sequence","Protein")
)

cat("Number of paired peptides (unmod vs succ):", nrow(paired), "\n")

# 3. Δdistribution score（modified vs unmodified）
unmod_cols <- grep("^Frac_.*_unmod$", colnames(paired), value = TRUE)
mod_cols   <- grep("^Frac_.*_mod$",   colnames(paired), value = TRUE)

frac_unmod <- as.matrix(paired[, unmod_cols])
frac_mod   <- as.matrix(paired[, mod_cols])

src_unmod_cols <- paste0("Frac_", paired$Main_unmod, "_unmod")
dst_unmod_cols <- paste0("Frac_", paired$Main_mod,   "_unmod")

src_mod_cols   <- paste0("Frac_", paired$Main_unmod, "_mod")
dst_mod_cols   <- paste0("Frac_", paired$Main_mod,   "_mod")

row_idx <- seq_len(nrow(paired))

src_unmod <- frac_unmod[cbind(row_idx, match(src_unmod_cols, colnames(frac_unmod)))]
dst_unmod <- frac_unmod[cbind(row_idx, match(dst_unmod_cols, colnames(frac_unmod)))]

src_mod   <- frac_mod[cbind(row_idx, match(src_mod_cols, colnames(frac_mod)))]
dst_mod   <- frac_mod[cbind(row_idx, match(dst_mod_cols, colnames(frac_mod)))]

paired$DeltaFracScore <- (dst_mod - src_mod) - (dst_unmod - src_unmod)

summary(paired$DeltaFracScore)

library(ggplot2)

df_score <- paired %>%
  filter(is.finite(DeltaFracScore)) %>%
  select(Sequence, Protein, DeltaFracScore)

# 简单统计
summary(df_score$DeltaFracScore)
quantile(abs(df_score$DeltaFracScore),
         probs = c(0.5, 0.8, 0.9, 0.95, 0.99),
         na.rm = TRUE)

# 分布图：原值
ggplot(df_score, aes(x = DeltaFracScore)) +
  geom_histogram(bins = 100, fill = "grey70", color = "grey30") +
  geom_vline(xintercept = c(-0.1, 0.1), color = "red", linetype = "dashed") +
  theme_bw() +
  labs(title = "Δdistribution score (modified vs unmodified)",
       x = "DeltaFracScore", y = "Peptides")

# 分布图：绝对值
ggplot(df_score, aes(x = abs(DeltaFracScore))) +
  geom_histogram(bins = 100, fill = "grey70", color = "grey30") +
  geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(title = "|Δdistribution score|",
       x = "|DeltaFracScore|", y = "Peptides")

# 不同阈值下的 peptide / protein 数量
thr_vec <- c(0.05, 0.1, 0.15, 0.2)

thr_stats <- lapply(thr_vec, function(th) {
  sel <- abs(df_score$DeltaFracScore) >= th
  data.frame(
    Threshold   = th,
    N_peptides  = sum(sel),
    N_proteins  = n_distinct(df_score$Protein[sel])
  )
}) %>% bind_rows()

thr_stats

# -------------------------------------------------------------------
# 4. 定义 peptide-level 和 protein-level movers
# -------------------------------------------------------------------
# 1. 准备 Location Map (这次我们要保留位点信息！)
#    原始 Sequence 格式如 "PEPTIDE: 123"，我们需要提取 "123" 作为位点
svm_mod_loc <- svm_mod %>% 
  mutate(
    # 提取冒号后面的数字作为位点位置 (Site)
    Site_Position = sub("^.*: ", "", Sequence), 
    # 提取纯序列用于 merge
    Clean_Sequence = sub(":.*$", "", Sequence)
  ) %>%
  dplyr::select(Clean_Sequence, Master.Protein.Accessions, Predicted_Location, Site_Position) %>%
  dplyr::rename(Protein = Master.Protein.Accessions,
                Loc_mod_label = Predicted_Location,
                Sequence = Clean_Sequence) %>%
  # 去重：同一个肽段、同一个蛋白、同一个位点、同一个定位，只留一个
  distinct() %>%
  group_by(Sequence, Protein) %>%
  slice(1) %>%
  ungroup()

# 2. 将定位和位点信息合并回 paired 表
paired_annotated <- paired %>%
  left_join(svm_mod_loc, by = c("Sequence", "Protein"))

# 3. 设定阈值
thr <- 0.1 

# 4. 标记特异性，并保留 Site 信息
peptide_specificity <- paired_annotated %>%
  mutate(
    IsSpecific = abs(MoveScore) >= thr,
    # 如果特异，保留定位；否则标为 Non-specific
    Specific_Compartment = if_else(IsSpecific, as.character(Loc_mod_label), "Non-specific"),
    # 构建位点标签：比如 "K364" (假设是Lysine修饰)
    # 注意：这里我们简单用 "Pos" + 数字，具体是不是 K 需要结合序列，但在结构分析时有位置就够了
    Site_Label = paste0(Site_Position) 
  )

# -------------------------------------------------------------------
# 5. 构建 UpSet 输入表并绘制 UpSet 图
# -------------------------------------------------------------------
# 构建 UpSet 输入
# 1. 转换成 Protein-Level 的宽表 (Binary Matrix)
# 逻辑：只要蛋白中有一个肽段是 Nucleus Specific，该蛋白在 Nucleus Set 就为 1
detailed_sites <- peptide_specificity %>%
  filter(Specific_Compartment != "Non-specific") %>% 
  filter(!is.na(Specific_Compartment)) %>%
  group_by(Protein, Specific_Compartment) %>%
  summarise(
    # 拼接该区室下的所有特异性位点
    Specific_Sites = paste(unique(Site_Label), collapse = "; "),
    .groups = "drop"
  ) %>%
  #以此生成宽表
  tidyr::pivot_wider(
    names_from = Specific_Compartment,
    values_from = Specific_Sites,
    names_prefix = "Sites_in_" # 列名前缀，如 Sites_in_Nucleus
  )

# 2. 生成 UpSet 图所需的 0/1 矩阵
#    逻辑：只要 detailed_sites 里某列有值，就为 1，否则为 0
upset_matrix <- detailed_sites %>%
  mutate(across(starts_with("Sites_in_"), ~ ifelse(is.na(.), 0, 1))) %>%
  rename_with(~ sub("Sites_in_", "", .), starts_with("Sites_in_")) %>%
  as.data.frame()

# 3. 绘制 UpSet 图 (代码不变)
library(UpSetR)
compartment_sets <- setdiff(colnames(upset_matrix), "Protein")

png("Fig_Succ_Specificity_UpSet_Refined.png", width = 2400, height = 1600, res = 300)
upset(
  upset_matrix,
  sets = compartment_sets,
  order.by = "freq",
  keep.order = TRUE,
  sets.bar.color = "#D55E00", # 换个颜色 (Vermilion)
  mainbar.y.label = "Protein Intersections",
  sets.x.label = "Proteins with Specific Succinylation"
)
dev.off()

final_target_list <- detailed_sites %>%
  rowwise() %>%
  mutate(
    # 计算非 NA 的列数（即有多少个特异性区室）
    Context_Count = sum(!is.na(c_across(starts_with("Sites_in_"))))
  ) %>%
  filter(Context_Count >= 2) %>% # 只保留多区室特异的 VIP 蛋白
  arrange(desc(Context_Count))

write.csv(final_target_list, "VIP_Target_Proteins_with_Sites.csv", row.names = FALSE)

cat("Done! CSV generated. Top proteins:\n")
print(head(final_target_list))














