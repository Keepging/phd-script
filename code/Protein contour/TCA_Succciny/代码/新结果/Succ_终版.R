library(MSnbase)
library(pRoloc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(data.table)

## ---------------------------------------------------------------------
## 0. 路径 & 输入文件
## ---------------------------------------------------------------------

setwd("D:\\博士\\TCA_Succciny\\代码")  # 若路径变了，这里改一下

mod_file   <- "succ_pep.rds"     # 修饰肽段 MSnSet
unmod_file <- "unmod_pep.rds"    # 未修饰肽段 MSnSet

## ---------------------------------------------------------------------
## 1. 加载原始数据
## ---------------------------------------------------------------------

cat(">>> 读取原始 MSnSet 对象...\n")
mod_msn   <- readRDS(mod_file)
unmod_msn <- readRDS(unmod_file)

## ---------------------------------------------------------------------
## 2. 未修饰组 SVM：确定蛋白“老家”
## ---------------------------------------------------------------------

cat(">>> 运行未修饰组 SVM，用于确定蛋白老家...\n")

# 简单过滤 Marker（与之前逻辑保持一致）
min_markers <- 20
unmod_counts <- table(fData(unmod_msn)$markers)
unmod_keep <- names(unmod_counts)[unmod_counts >= min_markers & names(unmod_counts) != "unknown"]
unmod_qc   <- unmod_msn[fData(unmod_msn)$markers %in% c(unmod_keep, "unknown"), ]

# SVM 参数沿用你之前的经验
unmod_model <- svmClassification(unmod_qc, cost = 16, sigma = 0.125)

# 提取结果
unmod_res <- data.frame(exprs(unmod_model))
unmod_res$Sequence <- rownames(unmod_res)
unmod_res$Predicted_Location <- fData(unmod_model)$svm
unmod_res$Score <- fData(unmod_model)$svm.scores
unmod_res$Master.Protein.Accessions <- fData(unmod_model)$Master.Protein.Accessions

# 高置信肽段
high_conf_unmod <- unmod_res %>% filter(Score > 0.7)

cat(">>> 未修饰组高置信肽段数：", nrow(high_conf_unmod), "\n")

## ---------------------------------------------------------------------
## 3. 修饰组 SVM：得到修饰后位置
## ---------------------------------------------------------------------

cat(">>> 运行修饰组 SVM...\n")

mod_counts <- table(fData(mod_msn)$markers)
mod_keep <- names(mod_counts)[mod_counts >= min_markers & names(mod_counts) != "unknown"]
mod_qc   <- mod_msn[fData(mod_msn)$markers %in% c(mod_keep, "unknown"), ]

mod_model <- svmClassification(mod_qc, cost = 16, sigma = 0.1)

mod_res <- data.frame(exprs(mod_model))
mod_res$Sequence <- rownames(mod_res)
mod_res$Predicted_Location <- fData(mod_model)$svm
mod_res$Score <- fData(mod_model)$svm.scores
mod_res$Master.Protein.Accessions <- fData(mod_model)$Master.Protein.Accessions

high_conf_mod <- mod_res %>% filter(Score > 0.7)

cat(">>> 修饰组高置信肽段数：", nrow(high_conf_mod), "\n")

## ---------------------------------------------------------------------
## 4. 构建“老家”映射：unmod_home_map
## ---------------------------------------------------------------------

cat(">>> 构建蛋白老家 (Unmod_Home_Loc)...\n")

unmod_home_map <- high_conf_unmod %>%
  group_by(Master.Protein.Accessions) %>%
  summarise(
    Unmod_Home_Loc = names(which.max(table(Predicted_Location))),
    Unmod_Pep_Count = n(),
    .groups = "drop"
  )

cat(">>> 老家映射的蛋白数：", nrow(unmod_home_map), "\n")

## ---------------------------------------------------------------------
## 5. 肽段层级：标记修饰肽段是否“移动”
##    - 只考虑：能找到老家的蛋白
## ---------------------------------------------------------------------

cat(">>> 标记肽段层面是否移动...\n")

peptide_level <- high_conf_mod %>%
  inner_join(unmod_home_map, by = "Master.Protein.Accessions") %>%
  mutate(
    Is_Moving_Peptide = (Predicted_Location != Unmod_Home_Loc),
    # 增加方向列：Format as "Home->New"
    Direction = paste0(Unmod_Home_Loc, "->", Predicted_Location)
  ) %>%
  dplyr::select(Sequence, Master.Protein.Accessions, Unmod_Home_Loc, 
         Mod_Location = Predicted_Location, Mod_Score = Score, 
         Is_Moving_Peptide, Direction)

## ---------------------------------------------------------------------
## 6. 蛋白层级：标记蛋白是否有“至少一个移动肽段”
## ---------------------------------------------------------------------

cat(">>> 标记蛋白层面是否移动...\n")

protein_level <- peptide_level %>%
  group_by(Master.Protein.Accessions, Unmod_Home_Loc) %>%
  summarise(
    Any_Moving_Peptide = any(Is_Moving_Peptide),
    .groups = "drop"
  )

## ---------------------------------------------------------------------
## 7. Flow 表：按蛋白汇总老家→修饰后位置
##    - 只统计真正“移动”的蛋白
## ---------------------------------------------------------------------

cat(">>> 构建蛋白层面的 Flow 表...\n")

flow_protein <- peptide_level %>%
  filter(Is_Moving_Peptide) %>%
  group_by(Master.Protein.Accessions, Unmod_Home_Loc, Mod_Location) %>%
  summarise(
    Max_Mod_Score = max(Mod_Score),
    .groups = "drop"
  )

## ---------------------------------------------------------------------
## 8. 输出核心结果表（供查阅/后续分析）
## ---------------------------------------------------------------------

cat(">>> 正在生成任务一要求的‘琥珀酰化移动肽段清洗表’...\n")

# 筛选条件：必须是移动的 (TRUE)
clean_moving_succi <- peptide_level %>%
  filter(Is_Moving_Peptide == TRUE)

# 导出 CSV
write.csv(clean_moving_succi, "Task1_Succi_Strict_Movers_Detailed.csv", row.names = FALSE)

cat(sprintf(">>> 完成。共筛选出 %d 个发生明确移动的琥珀酰化肽段。\n", nrow(clean_moving_succi)))
cat(">>> 文件已保存为: Task1_Succi_Strict_Movers_Detailed.csv\n")

cat(">>> 导出核心结果表...\n")

write.csv(peptide_level, "Succ_Peptide_Level_With_MoveFlag.csv", row.names = FALSE)
write.csv(protein_level, "Succ_Protein_Level_With_MoveFlag.csv", row.names = FALSE)
write.csv(flow_protein, "Succ_Protein_Flow_Table.csv", row.names = FALSE)

## =====================================================================
## 9. 绘图部分
## =====================================================================

## 9.1 肽段层面柱状图：每个区室移动肽段 / 总肽段
cat(">>> 绘制肽段层面柱状图...\n")

pep_stats_wide <- peptide_level %>%
  group_by(Unmod_Home_Loc) %>%
  summarise(
    Total_Peptides  = n(),
    Moving_Peptides = sum(Is_Moving_Peptide),
    .groups = "drop"
  ) %>%
  mutate(
    Stable_Peptides = Total_Peptides - Moving_Peptides,
    Move_Fraction   = Moving_Peptides / Total_Peptides,
    Stable_Fraction = Stable_Peptides / Total_Peptides,
    Label           = paste0(Moving_Peptides, " / ", Total_Peptides)
  )

pep_stats_long <- pep_stats_wide %>%
  select(Unmod_Home_Loc, Moving_Peptides, Stable_Peptides, 
         Move_Fraction, Stable_Fraction) %>%
  pivot_longer(
    cols = c(Moving_Peptides, Stable_Peptides, 
             Move_Fraction, Stable_Fraction),
    names_to  = c("Status", ".value"),
    names_sep = "_"
  )

p_pep_bar <- ggplot(pep_stats_long,
                    aes(x = Unmod_Home_Loc, y = Fraction, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = pep_stats_wide,
            aes(x = Unmod_Home_Loc, y = 1.02, label = Label),
            inherit.aes = FALSE, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  scale_fill_manual(
    values = c("Moving" = "grey40", "Stable" = "grey80"),
    labels = c("Moving" = "Moving Peptides",
               "Stable" = "Stable Peptides")
  )+
  labs(
    title    = "Peptide-Level Movement per Compartment",
    subtitle = "Proportion of succinylated peptides that change predicted location",
    x        = "Original Compartment (by Unmodified Protein Home)",
    y        = "Proportion of Peptides"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_pep_bar)
ggsave("Fig_Succ_Peptide_Movement_Proportion.png", p_pep_bar, width = 8, height = 5)





## 9.2 蛋白层面柱状图：每个“老家”总蛋白 / 移动蛋白

cat(">>> 绘制蛋白层面柱状图...\n")

# wide 版：一行一个 compartment，方便算比例和写标签
prot_stats_wide <- protein_level %>%
  group_by(Unmod_Home_Loc) %>%
  summarise(
    Total_Proteins  = n(),
    Moving_Proteins = sum(Any_Moving_Peptide),
    .groups = "drop"
  ) %>%
  mutate(
    Stable_Proteins = Total_Proteins - Moving_Proteins,
    Move_Fraction   = Moving_Proteins / Total_Proteins,
    Stable_Fraction = Stable_Proteins / Total_Proteins,
    Label           = paste0(Moving_Proteins, " / ", Total_Proteins)
  )

# long 版：专门用来画 stacked bar（y 已经是比例）
prot_stats_long <- prot_stats_wide %>%
  select(Unmod_Home_Loc, Moving_Proteins, Stable_Proteins, 
         Move_Fraction, Stable_Fraction) %>%
  pivot_longer(
    cols = c(Moving_Proteins, Stable_Proteins, 
             Move_Fraction, Stable_Fraction),
    names_to  = c("Status", ".value"),
    names_sep = "_"
  )

p_prot_bar <- ggplot(prot_stats_long,
                     aes(x = Unmod_Home_Loc, y = Fraction, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  # 在柱子上方标“Moving / Total”
  geom_text(data = prot_stats_wide,
            aes(x = Unmod_Home_Loc, y = 1.02, label = Label),
            inherit.aes = FALSE, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
  scale_fill_manual(
    values = c("Moving" = "grey40", "Stable" = "grey70"),
    labels = c("Moving" = "Moving Proteins",
               "Stable" = "Stable Proteins")
  ) +
  labs(
    title    = " ",
    subtitle = "Proportion of home proteins that change predicted location",
    x        = "Original Compartment (Unmodified Home)",
    y        = "Proportion of Proteins"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

print(p_prot_bar)

ggsave("Fig_Succ_Protein_Movement_Proportion.png", p_prot_bar, width = 8, height = 5)




## 9.3 蛋白层面桑基图：老家 → 修饰后位置（distinct proteins）
cat(">>> 绘制蛋白层面桑基图...\n")

flow_counts <- flow_protein %>%
  group_by(Unmod_Home_Loc, Mod_Location) %>%
  summarise(Protein_Count = n_distinct(Master.Protein.Accessions), .groups = "drop")

p_sankey <- ggplot(flow_counts,
                   aes(axis1 = Unmod_Home_Loc, axis2 = Mod_Location, y = Protein_Count)) +
  geom_alluvium(aes(fill = Unmod_Home_Loc), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "grey50") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Unmodified Home", "Succinylated Location"), expand = c(.05, .05)) +
  labs(
    title = "Protein Translocation Landscape upon Succinylation",
    subtitle = "Flows from unmodified home compartments to succinylated locations",
    y = "Number of Proteins",
    x = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(p_sankey)
ggsave("Fig_Succ_Protein_Sankey.png", p_sankey, width = 10, height = 6)


## 9.4 Flow 条形图：Top N 老家→新位置的路径（数字版）
cat(">>> 绘制 Flow 条形图 (Top flows)...\n")

topN <- 8

top_flows <- flow_counts %>%
  arrange(desc(Protein_Count)) %>%
  mutate(
    Flow_Label = paste(Unmod_Home_Loc, "→", Mod_Location)
  ) %>%
  head(topN)

p_flow_bar <- ggplot(top_flows, aes(x = reorder(Flow_Label, Protein_Count), y = Protein_Count)) +
  geom_col(fill = "#3288BD") +
  # 新增：在每个柱子上标出具体蛋白数
  geom_text(aes(label = Protein_Count),
            hjust = -0.1, size = 3) +
  coord_flip(clip = "off") +
  expand_limits(y = max(top_flows$Protein_Count) * 1.1) +
  labs(
    title = "Top Protein Flows upon Succinylation",
    subtitle = "Number of proteins moving from home to new locations",
    x = "Flow (Home \u2192 New Location)",
    y = "Number of Proteins"
  ) +
  theme_bw() +
  theme(plot.margin = margin(5.5, 30, 5.5, 5.5))  # 右边留点空位给数字


print(p_flow_bar)
ggsave("Fig_Succ_Protein_Flow_TopN.png", p_flow_bar, width = 8, height = 5)

cat("\n>>> 琥珀酰 pipeline 已完成：SVM + 肽段/蛋白统计 + 4 张图 + 3 个结果表。\n")
