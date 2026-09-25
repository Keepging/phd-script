library(dplyr)
library(readr)

# ==============================================================================
# 1. 读取数据 (这一步必须运行成功，变量才会存在！)
# ==============================================================================
cat(">>> 读取数据中...\n")

# 请确保路径正确，如果报错“No such file”，请检查文件名或路径
succi_movers <- read.csv("D:\\博士\\TCA_Succciny\\代码\\Task1_Succi_Strict_Movers_Detailed.csv")
phos_movers  <- read.csv("D:\\博士\\TCA_Succciny\\代码\\Strict_Movers_List.csv")
id_map       <- read.csv("D:\\博士\\TCA_Succciny\\代码\\p_g.csv")

cat(">>> 数据读取完成！\n")

# ==============================================================================
# 2. ID 转换 (Gene -> Uniprot)
# ==============================================================================
cat(">>> [Step 2] Mapping Gene Symbols to Uniprot IDs...\n")

# 整理字典 (确保无重复)
id_map_clean <- id_map %>%
  dplyr::select(Uniprot_ID = p, Gene = gene) %>%
  distinct(Gene, .keep_all = TRUE)

# 给磷酸化数据加上 Uniprot ID
phos_ready <- phos_movers %>%
  inner_join(id_map_clean, by = "Gene") %>%
  dplyr::select(Uniprot_ID, 
                Gene,
                Phos_ID = ID, 
                Phos_Dir = Direction, 
                Phos_Timepoint = Timepoint,
                Phos_Score = Movement_Score)

cat(sprintf(">>> ID Mapping complete. %d Phos events successfully mapped to Uniprot.\n", nrow(phos_ready)))

# ==============================================================================
# 3. Crosstalk 匹配
# ==============================================================================
cat(">>> [Step 3] Identifying Crosstalk (Same Protein, Both Moving)...\n")

# 整理琥珀酰化数据 (这里我也加上了 dplyr:: 以防报错)
succi_ready <- succi_movers %>%
  dplyr::select(Uniprot_ID = Master.Protein.Accessions, 
                Succi_Seq = Sequence, 
                Succi_Dir = Direction,
                Succi_Home = Unmod_Home_Loc,
                Succi_Dest = Mod_Location)

# 核心：Inner Join
crosstalk_events <- inner_join(succi_ready, phos_ready, by = "Uniprot_ID", relationship = "many-to-many")

# ==============================================================================
# 4. 结果导出
# ==============================================================================
if (nrow(crosstalk_events) > 0) {
  # 排序让结果更好看
  crosstalk_events <- crosstalk_events %>%
    arrange(Uniprot_ID, Gene)
  
  write.csv(crosstalk_events, "D:\\博士\\TCA_Succciny\\代码\\Task2_Crosstalk_Results_v2.csv", row.names = FALSE)
  
  cat("\n=======================================================\n")
  cat(sprintf(">>> SUCCESSS! Found %d Crosstalk events involving %d proteins.\n", 
              nrow(crosstalk_events), length(unique(crosstalk_events$Uniprot_ID))))
  cat(">>> Results saved to: 'Task2_Crosstalk_Results_v2.csv'\n")
  cat("=======================================================\n")
  
  # 打印前几行看看
  print(head(crosstalk_events %>% dplyr::select(Gene, Uniprot_ID, Succi_Dir, Phos_Dir)))
  
} else {
  cat("\n>>> Warning: No overlap found. Please check if 'p_g.csv' covers the relevant proteins.\n")
}