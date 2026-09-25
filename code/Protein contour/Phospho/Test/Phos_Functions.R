# ==============================================================================
# 文件名: 00_Phos_Functions.R
# 功能: 存放所有核心算法、绘图函数和 API 工具
# ==============================================================================

# --- 0. 加载所有依赖包 ---
library(tidyverse)
library(readxl)
library(igraph)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggalluvial)
library(UpSetR)
library(httr)      # API爬虫用
library(jsonlite)  # API数据解析用

# --- 1. 位置计算核心函数 ---

# 1.1 计算主 Fraction (FR1-FR6)
get_main_location <- function(df, time_pattern) {
  sub_df <- df %>% dplyr::select(contains(time_pattern))
  if (ncol(sub_df) != 24) return(NULL)
  
  mat_data <- 2^(as.matrix(sub_df))   # 反 log2
  group_factor <- rep(paste0("FR", 1:6), each = 4)
  
  row_means <- t(apply(mat_data, 1, function(x) {
    tapply(x, group_factor, mean, na.rm = TRUE)
  }))
  
  results <- data.frame(
    ID = df[[1]],
    Main_Loc = apply(row_means, 1, function(x) {
      if (sum(x, na.rm = TRUE) == 0) return(NA)
      names(which.max(x))
    }),
    stringsAsFactors = FALSE
  )
  return(results)
}

# 1.2 简化位置为三大区室
simplify_loc <- function(x) {
  case_when(
    grepl("FR1|FR2", x) ~ "Cytoplasm",
    grepl("FR3|FR4", x) ~ "Membrane",
    grepl("FR5|FR6", x) ~ "Nucleus",
    TRUE ~ "Other"
  )
}

# 1.3 计算 Movement Score 所需的百分比
get_super_compartment_percent <- function(df, time_pattern, id_col = "Phos_Key") {
  sub_df <- df %>% dplyr::select(all_of(id_col), contains(time_pattern))
  if (ncol(sub_df) < 2) return(NULL)
  
  mat_log2 <- as.matrix(sub_df[ , -1, drop = FALSE])
  mat_lin  <- 2 ^ mat_log2
  col_names <- colnames(sub_df)[-1]
  
  frac_labels <- ifelse(grepl("FR1", col_names, ignore.case = TRUE), "FR1",
                        ifelse(grepl("FR2", col_names, ignore.case = TRUE), "FR2",
                               ifelse(grepl("FR3", col_names, ignore.case = TRUE), "FR3",
                                      ifelse(grepl("FR4", col_names, ignore.case = TRUE), "FR4",
                                             ifelse(grepl("FR5", col_names, ignore.case = TRUE), "FR5",
                                                    ifelse(grepl("FR6", col_names, ignore.case = TRUE), "FR6", NA))))))
  
  if (any(is.na(frac_labels))) stop("列名匹配失败，请检查 FR1-FR6 标签。")
  
  fr_means <- t(apply(mat_lin, 1, function(x) tapply(x, frac_labels, mean, na.rm = TRUE)))
  fr_means <- as.data.frame(fr_means)
  rownames(fr_means) <- sub_df[[id_col]]
  
  cyto    <- fr_means$FR1 + fr_means$FR2
  memb    <- fr_means$FR3 + fr_means$FR4
  nucleus <- fr_means$FR5 + fr_means$FR6
  total   <- cyto + memb + nucleus
  total[total == 0 | is.na(total)] <- NA
  
  data.frame(
    ID    = rownames(fr_means),
    Cyto  = cyto / total,
    Memb  = memb / total,
    Nuc   = nucleus / total,
    Gene  = sub("_.*", "", rownames(fr_means)),
    stringsAsFactors = FALSE
  )
}

# 1.4 计算 Movement Score (核心逻辑)
calc_movement_score_pair <- function(df_pct, tp_stim, tp_ctrl = "EGF_CTRL") {
  cyto_ctrl <- paste0("Cyto_", tp_ctrl); memb_ctrl <- paste0("Memb_", tp_ctrl); nuc_ctrl  <- paste0("Nuc_",  tp_ctrl)
  cyto_stim <- paste0("Cyto_", tp_stim); memb_stim <- paste0("Memb_", tp_stim); nuc_stim  <- paste0("Nuc_",  tp_stim)
  
  df <- df_pct %>%
    filter(!is.na(.data[[cyto_ctrl]]) & !is.na(.data[[cyto_stim]]) &
             !is.na(.data[[memb_ctrl]]) & !is.na(.data[[memb_stim]]) &
             !is.na(.data[[nuc_ctrl]])  & !is.na(.data[[nuc_stim]]))
  
  if (nrow(df) == 0) return(NULL)
  
  d_cyto <- df[[cyto_stim]] - df[[cyto_ctrl]]
  d_memb <- df[[memb_stim]] - df[[memb_ctrl]]
  d_nuc  <- df[[nuc_stim]]  - df[[nuc_ctrl]]
  
  mv_score <- abs(d_cyto) + abs(d_memb) + abs(d_nuc)
  
  direction <- character(length(mv_score))
  for (i in seq_along(mv_score)) {
    dc <- d_cyto[i]; dm <- d_memb[i]; dn <- d_nuc[i]
    inc <- c(Cyto = dc, Memb = dm, Nuc = dn)
    dec <- inc
    inc[inc < 0] <- NA; dec[dec > 0] <- NA
    
    if (all(is.na(inc)) || all(is.na(dec))) {
      direction[i] <- "Complex/Global"
    } else {
      to   <- names(which.max(inc))
      from <- names(which.min(dec))
      # === 修改点：使用 -> 代替 → 解决乱码 ===
      direction[i] <- paste0(from, "->", to) 
    }
  }
  
  df %>% transmute(ID, Gene, Timepoint = tp_stim, 
                   Movement_Score = mv_score, Direction = direction,
                   d_Cyto = d_cyto, d_Memb = d_memb, d_Nuc = d_nuc)
}

# --- 2. 统计学计算函数 (p-value / q-value) ---

get_fr_reps_matrix <- function(df, tp, fr, id_col = "Phos_Key") {
  pattern <- paste0("^", tp, "_", fr, "_Rep[0-9]+$")
  cols <- grep(pattern, colnames(df), value = TRUE)
  if (length(cols) == 0) return(NULL)
  mat <- as.matrix(df[, cols, drop = FALSE])
  rownames(mat) <- df[[id_col]]
  return(mat)
}

get_compartment_p <- function(df, tp_stim, comp_name) {
  frs <- switch(comp_name, Cyto = c("FR1", "FR2"), Memb = c("FR3", "FR4"), Nuc  = c("FR5", "FR6"))
  mats_ctrl <- lapply(frs, function(fr) get_fr_reps_matrix(df, "EGF_CTRL", fr))
  mats_stim <- lapply(frs, function(fr) get_fr_reps_matrix(df, tp_stim,    fr))
  
  if (any(vapply(mats_ctrl, is.null, logical(1))) || any(vapply(mats_stim, is.null, logical(1)))) return(NULL)
  
  mat_ctrl <- Reduce(`+`, mats_ctrl)
  mat_stim <- Reduce(`+`, mats_stim)
  site_ids <- rownames(mat_ctrl)
  p_vec <- numeric(length(site_ids))
  
  for (i in seq_along(site_ids)) {
    x <- mat_ctrl[i, ]; y <- mat_stim[i, ]
    if (all(is.na(x)) || all(is.na(y))) { p_vec[i] <- NA_real_ } 
    else {
      tt <- try(t.test(x, y), silent = TRUE)
      p_vec[i] <- if (inherits(tt, "try-error")) NA_real_ else tt$p.value
    }
  }
  tibble(ID = site_ids, !!paste0("p_", comp_name, "_", tp_stim) := p_vec)
}

combine_p_fisher <- function(p_mat) {
  apply(as.matrix(p_mat), 1, function(p_row) {
    p_use <- p_row[!is.na(p_row)]
    if (length(p_use) == 0) return(NA_real_)
    stat <- -2 * sum(log(p_use))
    pchisq(stat, df = 2 * length(p_use), lower.tail = FALSE)
  })
}

# --- 3. Scope3P API 工具函数 ---
parse_id_info <- function(id_string) {
  parts <- strsplit(id_string, "_")[[1]]
  acc <- parts[1]
  pos_str <- gsub("[^0-9]", "", parts[2])
  return(list(acc = acc, pos = as.numeric(pos_str)))
}

fetch_scope3p_data <- function(acc, pos) {
  base_url <- "https://iomics.ugent.be/scop3p/api/modifications?accession="
  tryCatch({
    # 设置超时时间为 5 秒，防止卡死
    res <- GET(url = base_url, query = list(accession = acc), timeout(5)) 
    if (status_code(res) != 200) return("API_Error_Status")
    content_text <- content(res, "text", encoding = "UTF-8")
    json_data <- fromJSON(content_text)
    
    if (length(json_data) == 0 || nrow(json_data) == 0) return("Novel")
    match <- json_data %>% filter(position == pos)
    if (nrow(match) > 0) {
      func <- if("function" %in% colnames(match)) match[['function']][1] else "Known_No_Func"
      if (is.na(func) || func == "") func <- "Known_No_Func"
      return(func)
    } else {
      return("Novel") 
    }
  }, error = function(e) "Net_Fail")
}
# --- 4. 绘图与分析核心函数 (请追加到 00 号文件末尾) ---

analyze_dir_pair <- function(df, dir_pair, pair_name) {
  cat(sprintf("\n========== Processing: %s (%s <-> %s) ==========\n", 
              pair_name, dir_pair[1], dir_pair[2]))
  
  # 1. 从已有的严格名单中，筛选特定方向
  mv_pair <- df %>% filter(Direction %in% dir_pair)
  
  cat("Events in this pair:", nrow(mv_pair), "\n")
  if (nrow(mv_pair) == 0) {
    cat(">>> Skip: No movers found for this pair.\n")
    return(NULL)
  }
  
  # 2. 绘制 Barplot (时间分布)
  p_bar <- ggplot(mv_pair, aes(x = Timepoint, fill = Direction)) +
    geom_bar(position = "dodge") +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, position = position_dodge(0.9)) +
    labs(title = paste0(pair_name, " Movers Distribution"), y = "Count") +
    theme_bw()
  
  ggsave(sprintf("Fig_Strict_%s_Bar.png", pair_name), p_bar, width = 6, height = 5)
  
  # 3. 准备 UpSet 数据
  upset_df <- mv_pair %>%
    group_by(ID) %>%
    summarise(
      M_2min  = if(any(Timepoint == "EGF_2min")) 1 else 0,
      M_8min  = if(any(Timepoint == "EGF_8min")) 1 else 0,
      M_20min = if(any(Timepoint == "EGF_20min")) 1 else 0,
      M_90min = if(any(Timepoint == "EGF_90min")) 1 else 0
    ) %>% as.data.frame()
  
  # 4. 绘制 UpSet 图
  # 只在有数据时画图
  if(sum(upset_df[,-1]) > 0) {
    png_file <- sprintf("Fig_Strict_%s_UpSet.png", pair_name)
    png(png_file, width = 2000, height = 1500, res = 300)
    print(UpSetR::upset(upset_df, nsets = 4, order.by = "freq", 
                        mainbar.y.label = "Intersection Size", sets.x.label = "Set Size"))
    dev.off()
    cat(">>> Saved UpSet plot:", png_file, "\n")
  }
  
  # 5. 返回这一组的数据方便查看
  return(mv_pair)
}