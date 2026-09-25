library(dplyr)
library(readr)
library(httr)      
library(jsonlite)  

# ==============================================================================
# 0. 定义核心 API 函数 (语法修正版)
# ==============================================================================
fetch_scope3p_data <- function(acc, pos) {
  base_url <- "https://iomics.ugent.be/scop3p/api/modifications"
  
  tryCatch({
    res <- GET(
      url = base_url, 
      query = list(accession = acc), 
      add_headers("User-Agent" = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"),
      config(ssl_verifypeer = FALSE),
      timeout(30) # 给足 30 秒
    )
    
    if (status_code(res) != 200) {
      return(paste0("Error_", status_code(res))) 
    }
    
    # === 核心修改点：显式指定 as = "text" ===
    content_text <- content(res, as = "text", encoding = "UTF-8")
    
    json_data <- fromJSON(content_text)
    
    # 检查是否为空结果
    if (length(json_data) == 0) return("Novel")
    
    # 注意：有时候 json_data 是列表，有时候是 data.frame，这里做容错处理
    if (is.data.frame(json_data)) {
      if (nrow(json_data) == 0) return("Novel")
      df <- json_data
    } else if ("modifications" %in% names(json_data)) {
      # 处理你刚才贴出来的那种包含 modifications 字段的结构
      df <- as.data.frame(json_data$modifications)
    } else {
      return("Structure_Unknown")
    }
    
    # 查找特定位置
    # 确保 position 列存在
    if (!"position" %in% colnames(df)) return("Novel")
    
    match <- df %>% dplyr::filter(position == pos)
    
    if (nrow(match) > 0) {
      # 尝试获取功能描述，如果没有则返回 Known
      func <- if("function" %in% colnames(match)) match[['function']][1] else "Known_No_Func"
      if (is.na(func) || func == "") func <- "Known_No_Func"
      return(func)
    } else {
      return("Novel") 
    }
  }, error = function(e) {
    return(paste0("Fail: ", conditionMessage(e)))
  })
}

# ==============================================================================
# 1. 准备查询列表
# ==============================================================================
cat(">>> [Step 1] Preparing Scope3P Query List...\n")

phos_movers <- read.csv("D:\\博士\\TCA_Succciny\\代码\\Strict_Movers_List.csv")
id_map      <- read.csv("D:\\博士\\TCA_Succciny\\代码\\p_g.csv")

query_list <- phos_movers %>%
  dplyr::inner_join(id_map, by = c("Gene" = "gene")) %>%
  dplyr::rename(Uniprot_ID = p) %>%
  dplyr::mutate(Position = as.numeric(gsub(".*_[S|T|Y]([0-9]+)_.*", "\\1", ID))) %>%
  dplyr::filter(!is.na(Uniprot_ID), !is.na(Position)) %>%
  dplyr::distinct(Uniprot_ID, Position, .keep_all = TRUE)

total_sites <- nrow(query_list)
cat(sprintf(">>> Ready to query %d unique sites.\n", total_sites))

# ==============================================================================
# 2. 单个测试
# ==============================================================================
cat(">>> 测试连接中 (尝试第一个位点)...\n")
test_res <- fetch_scope3p_data(query_list$Uniprot_ID[1], query_list$Position[1])
cat(sprintf(">>> 测试结果: %s\n", test_res))

if (grepl("Fail", test_res)) {
  # 如果这次还失败，那就真的是见鬼了
  cat(">>> 依然报错，但网络看似是通的。请检查 httr 包版本。\n")
} else {
  cat(">>> 连接成功！开始批量运行...\n")
  
  # ==============================================================================
  # 3. 批量查询
  # ==============================================================================
  output_file <- "D:\\博士\\TCA_Succciny\\代码\\Task3_Scope3P_Results_v2.csv"
  if (file.exists(output_file)) file.remove(output_file) 
  
  for (i in 1:total_sites) {
    row <- query_list[i, ]
    status <- fetch_scope3p_data(row$Uniprot_ID, row$Position)
    
    res_row <- data.frame(
      Gene = row$Gene,
      Uniprot_ID = row$Uniprot_ID,
      Position = row$Position,
      Original_ID = row$ID,
      Scope3P_Status = status,
      Timestamp = format(Sys.time(), "%H:%M:%S")
    )
    
    write.table(res_row, output_file, sep = ",", 
                col.names = !file.exists(output_file), 
                row.names = FALSE, append = TRUE)
    
    if (i %% 5 == 0 || i == total_sites) {
      cat(sprintf("\r>>> Progress: %d/%d | Current: %s -> %s      ", 
                  i, total_sites, row$Gene, status))
      flush.console()
    }
    # 稍微慢一点，安全第一
    Sys.sleep(0.2) 
  }
  
  cat("\n\n>>> Summary:\n")
  final_data <- read.csv(output_file)
  print(table(final_data$Scope3P_Status))
  cat(sprintf(">>> Done! Full results saved to '%s'.\n", output_file))
}