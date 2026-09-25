# ==============================================================================
# 文件名: 01_Phos_Main_Analysis.R
# 功能: 磷酸化主线分析 - 算分、统计、严格筛选与绘图
# ==============================================================================

# 1. 加载工具箱
source("D:\\博士\\Phospho\\Test\\Phos_Functions.R")

# 2. 读取数据
data_path <- "D:\\博士\\Phospho\\Full_data.xlsx" # 请确认路径
cat(">>> [Step 1] Loading Data...\n")
df_prot <- read_excel(data_path, sheet = "S5B-Log2 proc HeLa+EGF PROT")
df_phos <- read_excel(data_path, sheet = "S5D-Log2 proc HeLa+EGF PHOS")
colnames(df_prot)[1] <- "Gene"
colnames(df_phos)[1] <- "Phos_Key"

# --- 核心计算模块 ---

cat(">>> [Step 2] Calculating Movement Score & q-values...\n")
timepoints_mv <- c("EGF_CTRL", "EGF_2min", "EGF_8min", "EGF_20min", "EGF_90min")

# 2.1 计算百分比
pct_list <- list()
for (tp in timepoints_mv) {
  tmp <- get_super_compartment_percent(df_phos, tp, id_col = "Phos_Key")
  if (!is.null(tmp)) {
    colnames(tmp)[2:4] <- paste0(c("Cyto", "Memb", "Nuc"), "_", tp)
    pct_list[[tp]] <- tmp
  }
}
pct_merged <- purrr::reduce(pct_list, dplyr::full_join, by = c("ID", "Gene"))

# 2.2 计算 Movement Score
mv_list <- list()
for (tp in setdiff(timepoints_mv, "EGF_CTRL")) {
  mv_list[[tp]] <- calc_movement_score_pair(pct_merged, tp_stim = tp)
}
movement_df <- bind_rows(mv_list)

# 2.3 计算 q-value
p_list_all <- list()
for (tp in setdiff(timepoints_mv, "EGF_CTRL")) {
  p_cyto <- get_compartment_p(df_phos, tp_stim = tp, comp_name = "Cyto")
  p_memb <- get_compartment_p(df_phos, tp_stim = tp, comp_name = "Memb")
  p_nuc  <- get_compartment_p(df_phos, tp_stim = tp, comp_name = "Nuc")
  
  p_df <- list(p_cyto, p_memb, p_nuc) %>% 
    purrr::discard(is.null) %>% 
    purrr::reduce(dplyr::full_join, by = "ID")
  if (nrow(p_df) == 0) next
  
  p_cols <- grep(paste0("p_.*_", tp, "$"), colnames(p_df), value = TRUE)
  
  p_df <- p_df %>%
    mutate(p_combined = combine_p_fisher(dplyr::select(., all_of(p_cols))),
           q_value    = p.adjust(p_combined, method = "BH"),
           Timepoint  = tp)
  p_list_all[[tp]] <- p_df %>% dplyr::select(ID, Timepoint, q_value)
}

# 2.4 合并最终总表
p_stats_all <- bind_rows(p_list_all)
movement_df <- movement_df %>% left_join(p_stats_all, by = c("ID", "Timepoint"))

# 保存包含所有数据的大表
write.csv(movement_df, "Phos_MovementScore_All_Data.csv", row.names = FALSE)

# --- 关键筛选模块 (Strict Filtering) ---

cat(">>> [Step 3] Filtering Strict Movers (Score >= 0.1 & q < 0.05)...\n")
strict_movers <- movement_df %>%
  filter(!is.na(q_value)) %>%
  filter(Movement_Score >= 0.1, q_value < 0.05)

# 保存这份“严格”的名单，供下一步 Scope3P 使用
# 注意：只保存唯一的 ID 即可，或者保存完整信息方便查阅
write.csv(strict_movers, "Strict_Movers_List.csv", row.names = FALSE)
cat(sprintf(">>> Filtered: %d events found. Saved to 'Strict_Movers_List.csv'.\n", nrow(strict_movers)))

# --- 绘图模块 (Barplot / UpSet) ---
# 只有当你已经生成了 strict_movers (严格筛选名单) 后，才能运行下面这些

cat(">>> [Step 4] Visualizing Specific Pairs...\n")

# 确保数据存在
if(!exists("strict_movers")) {
  strict_movers <- read.csv("Strict_Movers_List.csv") # 如果之前跑过，重新读取防止报错
}

# 1. 分析 细胞质 <-> 细胞核 (CytoNuc)
# 这一步会自动画出 Fig_Strict_CytoNuc_Bar.png 和 UpSet 图
res_cn <- analyze_dir_pair(strict_movers, 
                           dir_pair = c("Cyto→Nuc", "Nuc→Cyto"), 
                           pair_name = "CytoNuc")

# 2. 分析 细胞质 <-> 细胞膜 (CytoMemb)
res_cm <- analyze_dir_pair(strict_movers, 
                           dir_pair = c("Cyto→Memb", "Memb→Cyto"), 
                           pair_name = "CytoMemb")

# 3. 分析 细胞膜 <-> 细胞核 (MembNuc)
res_mn <- analyze_dir_pair(strict_movers, 
                           dir_pair = c("Memb→Nuc", "Nuc→Memb"), 
                           pair_name = "MembNuc")

cat("\n>>> All plots generated. Check your folder.\n")