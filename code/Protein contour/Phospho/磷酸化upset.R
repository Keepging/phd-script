library(dplyr)
library(ggplot2)
library(UpSetR)
library(grid)

# ==============================================================================
# 定义分析函数
# ==============================================================================
analyze_dir_pair <- function(df, dir_pair, pair_name) {
  
  cat(sprintf("\n========== Processing: %s (%s <-> %s) ==========\n",
              pair_name, dir_pair[1], dir_pair[2]))
  
  # --- 1. 数据过滤与准备 ---
  # 确保数值列类型正确
  df$q_value <- as.numeric(as.character(df$q_value))
  df$Movement_Score <- as.numeric(as.character(df$Movement_Score))
  
  # 筛选符合条件的移动行
  mv_pair <- df %>%
    filter(Direction %in% dir_pair) %>%
    filter(!is.na(q_value), !is.na(Movement_Score)) %>% 
    filter(Movement_Score >= 0.1, q_value < 0.05)
  
  cat("Number of moving rows (filtered):", nrow(mv_pair), "\n")
  
  if (nrow(mv_pair) == 0) {
    cat(">>> STOP: No movers found. Skipping.\n")
    return(invisible(NULL))
  }
  
  # --- 2. 柱状图 (Barplot) ---
  dir_time <- mv_pair %>%
    count(Timepoint, Direction, name = "Count")
  
  p_dir_time <- ggplot(dir_time,
                       aes(x = Timepoint, y = Count, fill = Direction)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = Count),
              position = position_dodge(width = 0.9),
              vjust = -0.3, size = 3) +
    labs(title = sprintf("%s movers over time", pair_name)) +
    theme_bw()
  
  # 保存为 PNG (白色背景)
  bar_file <- sprintf("Fig_Phos_Movement_%s_Bar.png", pair_name)
  png(bar_file, width = 2000, height = 1600, res = 300, bg = "white")
  print(p_dir_time)
  dev.off()
  cat("Saved Barplot:", bar_file, "\n")
  
  # --- 3. UpSet 图 (核心修正部分) ---
  
  # 清理 Timepoint 字符串，防止空格干扰
  mv_pair$Timepoint <- trimws(as.character(mv_pair$Timepoint))
  
  # 生成 0/1 矩阵 (使用 if/else 确保无 NA)
  overlap_df <- mv_pair %>%
    group_by(ID) %>%
    summarise(
      Gene        = first(Gene),
      Moves_2min  = if(any(Timepoint == "EGF_2min")) 1 else 0,
      Moves_8min  = if(any(Timepoint == "EGF_8min")) 1 else 0,
      Moves_20min = if(any(Timepoint == "EGF_20min")) 1 else 0,
      Moves_90min = if(any(Timepoint == "EGF_90min")) 1 else 0,
      .groups     = "drop"
    )
  
  # 保存 CSV
  csv_overlap <- sprintf("Phos_Movement_%s_Overlap.csv", pair_name)
  write.csv(overlap_df, csv_overlap, row.names = FALSE)
  
  # 提取矩阵用于绘图
  upset_input <- overlap_df %>%
    select(Moves_2min, Moves_8min, Moves_20min, Moves_90min) %>%
    as.data.frame() # 必须转换为 data.frame
  
  # 强制转换为数值型 (防止 integer 导致某些版本不兼容)
  upset_input[] <- lapply(upset_input, as.numeric)
  
  # 诊断信息
  cat("--- UpSet Data Check ---\n")
  print(head(upset_input))
  cat("Column Sums:", colSums(upset_input), "\n")
  cat("Total Events:", sum(upset_input), "\n")
  
  if (sum(upset_input) == 0) {
    cat(">>> ERROR: Matrix is empty (all zeros). Cannot plot UpSet.\n")
  } else {
    
    # === 方法 A: 保存为 PDF (最稳妥，用于检查内容) ===
    # PDF 是矢量图，不受分辨率影响，如果这个能看，说明数据没问题
    pdf_file <- sprintf("Fig_Phos_Movement_%s_UpSet.pdf", pair_name)
    pdf(pdf_file, width = 8, height = 6) # 标准尺寸英寸
    print(
      UpSetR::upset(upset_input, 
                    nsets = 4, 
                    mainbar.y.label = "Intersection Size",
                    sets.x.label = "Set Size",
                    order.by = "freq")
    )
    dev.off()
    cat("Saved UpSet (PDF):", pdf_file, "<-- 请优先检查这个文件!\n")
    
    # === 方法 B: 保存为高分辨率 PNG (已调整字体大小) ===
    png_file <- sprintf("Fig_Phos_Movement_%s_UpSet.png", pair_name)
    png(png_file, width = 2400, height = 1800, res = 300, bg = "white")
    
    print(
      UpSetR::upset(upset_input, 
                    nsets = 4, 
                    order.by = "freq",
                    point.size = 3.5,    # 放大点
                    line.size = 1.2,     # 加粗线
                    text.scale = 1.5)    # 放大文字 (关键! 否则res=300时字太小)
    )
    dev.off()
    cat("Saved UpSet (PNG):", png_file, "\n")
  }
  
  # --- 4. 双向位点导出 ---
  bidir <- mv_pair %>%
    group_by(ID, Gene) %>%
    summarise(
      Has_dir1 = any(Direction == dir_pair[1]),
      Has_dir2 = any(Direction == dir_pair[2]),
      .groups  = "drop"
    ) %>%
    filter(Has_dir1 & Has_dir2)
  
  write.csv(bidir, sprintf("Phos_Movement_%s_Bidir.csv", pair_name), row.names = FALSE)
  
  invisible(list(upset_input = upset_input))
}

# ==============================================================================
# 运行分析 (请确保 environment 中有 movement_df)
# ==============================================================================

if (!exists("movement_df")) {
  stop("Error: 'movement_df' not found. Please load your data first.")
}

cat("Starting Analysis...\n")

# 运行第一组
res_CytoNuc  <- analyze_dir_pair(movement_df, c("Cyto→Nuc","Nuc→Cyto"), "CytoNuc")

# 如果需要运行其他组，取消下面的注释
# res_CytoMemb <- analyze_dir_pair(movement_df, c("Cyto→Memb","Memb→Cyto"), "CytoMemb")
# res_MembNuc  <- analyze_dir_pair(movement_df, c("Memb→Nuc","Nuc→Memb"), "MembNuc")

cat("\nDone.\n")


# ==============================================================================
# 运行另外两个部分的 TIFF 图生成
# ==============================================================================

# 确保函数已经定义（如果还没定义，请先运行上一次回复中的函数定义部分）

if (!exists("movement_df")) {
  stop("Error: 'movement_df' not found. Please load your data first.")
}

cat("\n>>> Starting Remaining Analyses (TIFF Output)...\n")

# --- 1. 处理 CytoMemb ---
res_CytoMemb <- analyze_dir_pair(movement_df, 
                                 dir_pair  = c("Cyto→Memb","Memb→Cyto"), 
                                 pair_name = "CytoMemb")

if (!is.null(res_CytoMemb)) {
  tiff_file <- "Fig_Phos_Movement_CytoMemb_UpSet.tiff"
  tiff(tiff_file, width = 2400, height = 1800, res = 300, compression = "lzw", bg = "white")
  
  # 显式使用 print 确保渲染
  print(UpSetR::upset(
    res_CytoMemb$upset_input, 
    nsets           = 4, 
    order.by        = "freq",
    point.size      = 3.5, 
    line.size       = 1.2, 
    text.scale      = c(2, 2, 1.5, 1.5, 2, 1.5), # 放大所有文字标签
    sets.x.label    = "CytoMemb sites",
    mainbar.y.label = "Intersections"
  ))
  
  dev.off()
  cat("Successfully saved TIFF:", tiff_file, "\n")
}

# --- 2. 处理 MembNuc ---
res_MembNuc <- analyze_dir_pair(movement_df, 
                                dir_pair  = c("Memb→Nuc","Nuc→Memb"), 
                                pair_name = "MembNuc")

if (!is.null(res_MembNuc)) {
  tiff_file <- "Fig_Phos_Movement_MembNuc_UpSet.tiff"
  tiff(tiff_file, width = 2400, height = 1800, res = 300, compression = "lzw", bg = "white")
  
  print(UpSetR::upset(
    res_MembNuc$upset_input, 
    nsets           = 4, 
    order.by        = "freq",
    point.size      = 3.5, 
    line.size       = 1.2, 
    text.scale      = c(2, 2, 1.5, 1.5, 2, 1.5), 
    sets.x.label    = "MembNuc sites",
    mainbar.y.label = "Intersections"
  ))
  
  dev.off()
  cat("Successfully saved TIFF:", tiff_file, "\n")
}

cat("\n>>> All processing done. Please check your folder for TIFF files.\n")