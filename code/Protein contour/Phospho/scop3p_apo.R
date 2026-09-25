library(tidyverse)
library(httr)
library(jsonlite)

# ==============================================================================
# 1. 最终修正的 API 查询函数 (基于 /api/modifications 接口)
# ==============================================================================
check_scop3p_final <- function(uniprot_acc, site_num) {
  
  # 关键修改：切换到 Image 1 中的 /api/modifications 接口
  # 这个接口直接返回蛋白层面的位点，而不是肽段层面，更准确且不易出错
  url <- paste0("https://iomics.ugent.be/scop3p/api/modifications?accession=", uniprot_acc)
  
  # 发送请求
  resp <- tryCatch(GET(url, timeout(10)), error = function(e) NULL)
  
  # 错误处理
  if (is.null(resp) || status_code(resp) != 200) return("API_Error")
  
  # 解析 JSON
  content_txt <- content(resp, "text", encoding = "UTF-8")
  if (nchar(content_txt) < 5) return("NoData") 
  
  json_dat <- tryCatch(fromJSON(content_txt), error = function(e) NULL)
  if (is.null(json_dat)) return("NoData")
  
  # --- 核心提取逻辑 ---
  # 根据图1，返回的 JSON 应该包含一个 "modifications" 字段
  # 这是一个数据框，里面有一列叫 "position"
  
  known_sites <- c()
  
  if (!is.null(json_dat$modifications)) {
    # 提取 position 列
    if ("position" %in% colnames(json_dat$modifications)) {
      known_sites <- json_dat$modifications$position
    }
  }
  
  # 确保转为数字进行比较
  if (as.numeric(site_num) %in% as.numeric(known_sites)) {
    return("Known")
  } else {
    return("Novel")
  }
}

# ==============================================================================
# 2. 现场验证 (P12270, Site 1187)
# ==============================================================================
cat(">>> 正在进行最终验证 (P12270, Site 1187)...\n")
# 应该返回 "Known"
test_result <- check_scop3p_final("P12270", "1187") 

cat(sprintf(">>> 测试结果: %s\n", test_result))

if(test_result == "Known") {
  cat("✅ 验证成功！API 接口已正确切换。\n")
} else {
  stop("❌ 验证失败。请检查网络或 Uniprot ID。")
}

# ==============================================================================
# 3. 全量运行 (如果 target_list 在内存中)
# ==============================================================================
# 如果你已经有 target_list，直接运行下面这段更新状态
if (exists("target_list")) {
  cat("\n>>> 开始全量更新数据库状态...\n")
  pb <- txtProgressBar(min = 0, max = nrow(target_list), style = 3)
  
  for (i in 1:nrow(target_list)) {
    # 只更新之前没有结果的，或者为了保险起见全部重跑
    # 这里我们全部重跑，确保准确
    if (!is.na(target_list$Uniprot_ID[i]) && !is.na(target_list$Site_Num[i])) {
      
      status <- check_scop3p_final(target_list$Uniprot_ID[i], target_list$Site_Num[i])
      target_list$Scop3P_Status[i] <- status
      
      Sys.sleep(0.05) # 稍微快一点也没事，这个接口响应应该很快
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 打印最终统计
  cat("\n>>> 最终 Scop3P 统计结果：\n")
  print(table(target_list$Scop3P_Status))
  
  # 保存
  write.csv(target_list, "Final_Result_Scop3P_Corrected_API.csv", row.names = FALSE)
  cat("\n>>> 结果已保存。\n")
}

library(tidyverse)
library(ggplot2)

# 假设 target_list 是你刚刚跑完的 712/154 那个表
# target_list <- read.csv("Final_Result_Scop3P_Corrected_API.csv")

# ==============================================================================
# 1. 最终饼图 (High Confidence Validation)
# ==============================================================================
stats_df <- target_list %>%
  filter(Scop3P_Status %in% c("Known", "Novel")) %>%
  count(Scop3P_Status) %>%
  mutate(prop = n / sum(n) * 100) %>%
  mutate(label = paste0(Scop3P_Status, "\n", n, " (", round(prop, 1), "%)"))

# 配色策略：Known用灰色(表示背景)，Novel用亮红色(表示发现)
p_pie <- ggplot(stats_df, aes(x = "", y = prop, fill = Scop3P_Status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Known" = "grey80", "Novel" = "#E41A1C")) + 
  theme_void() +
  labs(title = "Database Validation (Scop3P)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Fig_Final_Validwwation_Pie.png", p_pie, width = 4, height = 4)

# ==============================================================================
# 2. 最终的 Top 10 Novel 列表 (The Golden List)
# ==============================================================================
# 这些是真正的 gems
top_novel_final <- target_list %>%
  filter(Scop3P_Status == "Novel") %>%
  arrange(desc(Prob_Score)) %>% # 按置信度排序
  select(Gene, Phos_Key, Group_Shift, Prob_Score) %>%
  head(10)

# 保存表格供 PPT 使用
write.csv(top_novel_final, "Table_Final_Noveeel_Candidates.csv", row.names = FALSE)

print(p_pie)
print(top_novel_final)




library(tidyverse)
library(ggplot2)

# ==============================================================================
# 1. 加载数据 (修复 "找不到对象" 的报错)
# ==============================================================================
# 这里的文件名必须是你之前跑完 API 后保存的那个文件名
# 如果你的文件名不一样，请在这里修改
csv_file <- "Final_Result_Scop3P_Corrected_API.csv"

if (file.exists(csv_file)) {
  cat(">>> 正在读取保存的 Scop3P 验证结果...\n")
  target_list <- read.csv(csv_file)
  cat(sprintf(">>> 读取成功！包含 %d 行数据。\n", nrow(target_list)))
} else {
  stop("错误：找不到文件 'Final_Result_Scop3P_Corrected_API.csv'。请确认你是否已经跑完了 API 查询步骤并保存了文件。")
}

# ==============================================================================
# 2. 总体饼图 (PPT大字版)
# ==============================================================================
cat(">>> 正在生成饼图 (大字版)...\n")

stats_df <- target_list %>%
  filter(Scop3P_Status %in% c("Known", "Novel")) %>%
  count(Scop3P_Status) %>%
  mutate(prop = n / sum(n) * 100) %>%
  mutate(label = paste0(Scop3P_Status, "\n", n, "\n(", round(prop, 1), "%)"))

p_pie <- ggplot(stats_df, aes(x = "", y = prop, fill = Scop3P_Status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  
  # [修改] 标签字号改为 8 (非常大)，并加粗
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), 
            size = 8, fontface = "bold", lineheight = 1.2) +
  
  scale_fill_manual(values = c("Known" = "grey80", "Novel" = "#E41A1C")) + 
  theme_void() +
  labs(title = "Overall Discovery Rate") +
  theme(
    # [修改] 标题字号改为 24，加粗
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    legend.position = "none" # 饼图里已经写了字，不需要图例占地方了
  )

# 保存时稍微调大一点尺寸，避免字挤出去
ggsave("Fig_Final_Validation_Pie_PPT.png", p_pie, width = 6, height = 6, dpi = 300)
print(p_pie)

# ==============================================================================
# 3. 柱状图 (修复版：智能标签防止重叠)
# ==============================================================================
cat(">>> 正在生成修复版柱状图...\n")

# 重新计算并生成智能标签
breakdown_data <- target_list %>%
  filter(Scop3P_Status %in% c("Known", "Novel")) %>%
  group_by(Group_Shift, Scop3P_Status) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Group_Shift) %>%
  mutate(
    Total = sum(Count),
    Percent = Count / Total * 100,
    
    # --- [核心修改] 智能标签逻辑 ---
    # 如果占比小于 12% 或者 数量小于 20，格子太矮，只显示数字，不显示百分比
    # 否则显示两行：数字 + 百分比
    Label = ifelse(Percent < 12 | Count < 20, 
                   as.character(Count), 
                   paste0(Count, "\n(", round(Percent, 0), "%)"))
  )

p_breakdown <- ggplot(breakdown_data, aes(x = Group_Shift, y = Count, fill = Scop3P_Status)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") +
  
  # --- [核心修改] 调整字体参数 ---
  geom_text(aes(label = Label), 
            position = position_stack(vjust = 0.5), 
            size = 5,           # 从 6 改为 5 (依然很大，但不会挤)
            lineheight = 0.8,   # 行距设为 0.8 (让两行字靠得更近)
            fontface = "bold", 
            color = "white") +
  
  scale_fill_manual(values = c("Known" = "grey70", "Novel" = "#E41A1C")) +
  
  labs(
    title = "Discovery Rate by Path", 
    x = NULL, 
    y = "Number of Sites",
    fill = "Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.text.x = element_text(face = "bold", size = 16, angle = 20, hjust = 1, color = "black"),
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    axis.title.y = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14),
    legend.position = "top", 
    panel.grid.major.x = element_blank()
  )

# 保存
ggsave("Fig_Scop3P_Breakdown_Bar_Fixed.png", p_breakdown, width = 10, height = 7, dpi = 300)
print(p_breakdown)

cat("\n>>> 修复版柱状图已生成！重叠问题已解决。\n")