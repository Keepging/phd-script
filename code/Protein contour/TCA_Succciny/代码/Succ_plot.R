cat("\n=== 开始进行蛋白层面归位分析 ===\n")

# A. 确定未修饰蛋白的“老家” (Consensus Location)
# 逻辑：一个蛋白可能有10条肽段，6条在 Cyto，4条在 Nuc -> 老家 = Cyto
unmod_home_map <- high_conf_unmod %>%
  group_by(Master.Protein.Accessions) %>%
  summarise(
    Unmod_Home_Loc = names(which.max(table(Predicted_Location))),
    Unmod_Pep_Count = n()
  )

# B. 将“老家”信息贴到修饰肽段上
# 逻辑：这是蛋白X的修饰肽段，它现在在哪？蛋白X的老家在哪？
merged_analysis <- high_conf_results %>%
  inner_join(unmod_home_map, by = "Master.Protein.Accessions") %>%
  select(Sequence, Master.Protein.Accessions, Predicted_Location, Unmod_Home_Loc)

# C. 找出“离家出走”的 (Translocated)
translocated_df <- merged_analysis %>%
  filter(Predicted_Location != Unmod_Home_Loc)

# ==============================================================================
# 3. 关键统计结果输出 (给师兄看的数字)
# ==============================================================================

cat("\n#######################################################\n")
cat("               FINAL STATISTICS REPORT                 \n")
cat("#######################################################\n")

# 1. 肽段层面
num_trans_peps <- nrow(translocated_df)
cat(sprintf("1. [Peptide Level] 发生易位的修饰肽段数量: %d\n", num_trans_peps))

# 2. 蛋白层面 (这就是那个 399/400 的来源验证)
num_trans_prots <- length(unique(translocated_df$Master.Protein.Accessions))
cat(sprintf("2. [Protein Level] 发生易位的蛋白数量 (Unique Proteins): %d\n", num_trans_prots))

cat("-------------------------------------------------------\n")
cat("结论话术建议：\n")
cat(sprintf("师兄，我重新核算了。我们发现了 %d 个修饰肽段发生了定位改变，\n", num_trans_peps))
cat(sprintf("这些肽段归属于 %d 个独立的蛋白。\n", num_trans_prots))
cat(sprintf("所以之前 PPT 里的数字应该是 %d (请填入上方Protein Level的数字)。\n", num_trans_prots))
cat("#######################################################\n")

# 保存结果表，防止师兄查岗
write.csv(translocated_df, "Final_Translocation_Candidates_ProteinLevel.csv", row.names = FALSE)


# ==============================================================================
# 4. 任务一修正：基于蛋白层面的区室重叠统计 (用于画图)
# ==============================================================================
# 之前的逻辑是查 clean_seq (全0)，现在我们查 Protein ID

cat("\n>>> 正在生成修正后的统计表 (Protein-Level Overlap)...\n")

compartments <- unique(c(high_conf_results$Predicted_Location, high_conf_unmod$Predicted_Location))
compartments <- compartments[!compartments %in% c("unknown", NA)]

stats_list <- list()

for (comp in compartments) {
  # 1. 该区室有哪些“常驻居民”蛋白 (基于未修饰数据)
  home_proteins <- unmod_home_map %>% 
    filter(Unmod_Home_Loc == comp) %>% 
    pull(Master.Protein.Accessions)
  
  count_unmod_prots <- length(home_proteins)
  
  # 2. 这些“常驻民”里，有多少人在“其他地方”被修饰了？
  # 在 merged_analysis 里找：蛋白是常驻民(home=comp) 且 修饰位置 != comp
  moved_proteins <- merged_analysis %>%
    filter(Master.Protein.Accessions %in% home_proteins) %>%
    filter(Predicted_Location != comp) %>%
    pull(Master.Protein.Accessions) %>%
    unique()
  
  count_moved <- length(moved_proteins)
  
  stats_list[[comp]] <- data.frame(
    Compartment = comp,
    Total_Unmod_Proteins = count_unmod_prots,
    Proteins_Found_Modified_Elsewhere = count_moved
  )
}

final_stats_df <- bind_rows(stats_list)
print(final_stats_df)

# ==============================================================================
# 5. 修正后的可视化 (修复了图例反了的问题)
# ==============================================================================

plot_df <- final_stats_df %>%
  mutate(Pure_Unmod = Total_Unmod_Proteins - Proteins_Found_Modified_Elsewhere) %>%
  select(Compartment, Proteins_Found_Modified_Elsewhere, Pure_Unmod) %>%
  pivot_longer(cols = c("Proteins_Found_Modified_Elsewhere", "Pure_Unmod"), 
               names_to = "Status", values_to = "Count")

# 修正图例逻辑：
# Pure_Unmod -> 灰色 (没有在别处发现修饰)
# Proteins_Found_Modified_Elsewhere -> 红色 (发现了异位修饰)

p_corrected <- ggplot(plot_df, aes(x = Compartment, y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c("Pure_Unmod" = "grey70", "Proteins_Found_Modified_Elsewhere" = "#D53E4F"),
    # 这里的 labels 顺序必须和 values 里的 key 对应，或者利用 breaks 指定
    breaks = c("Proteins_Found_Modified_Elsewhere", "Pure_Unmod"),
    labels = c("Modified Elsewhere (Protein Level)", "Always Unmod Here")
  ) +
  labs(title = "Protein-Level Compartment Specificity",
       subtitle = "Of proteins residing in X (Unmod), how many are modified elsewhere?",
       y = "Proportion of Proteins", x = "Subcellular Compartment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_corrected)
ggsave("Fig_Protein_Level_Overlap_Corrected.png", p_corrected, width = 8, height = 5)

cat("\n所有任务完成。请查看控制台输出的数字以及新生成的图片。\n")