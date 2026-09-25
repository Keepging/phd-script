# ==============================================================================
# 文件名: 03_DeepMVP_Prep.R
# 功能: 自动化准备 DeepMVP 网页版所需的批量上传文件
# 输入: Strict_Movers_List.csv
# 输出: DeepMVP_Upload_Sequences.fasta, DeepMVP_Batch_Mutation_Input.txt
# ==============================================================================

source("D:\\博士\\Phospho\\Test\\Phos_Functions.R") # 加载依赖
library(httr)
library(jsonlite) # 如果没有安装请运行 install.packages("jsonlite")

# --- 1. 读取严格筛选名单 ---
if(!file.exists("Strict_Movers_List.csv")) {
  stop("找不到 'Strict_Movers_List.csv'")
}
movers <- read.csv("Strict_Movers_List.csv")
unique_ids <- unique(movers$ID)
cat(sprintf(">>> 共有 %d 个唯一位点需要分析。\n", length(unique_ids)))

# --- 2. 解析并转换 ID (Gene Name -> Uniprot Accession) ---
# 定义一个函数，通过 API 将基因名转为 ID
get_accession_from_gene <- function(gene_name) {
  # 搜索 Swiss-Prot (reviewed=true) 且 人类 (9606)
  # 如果你的数据是小鼠，请把 9606 改为 10090
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  query <- paste0("gene_exact:", gene_name, " AND reviewed:true AND organism_id:9606")
  
  tryCatch({
    res <- GET(base_url, query = list(query = query, format = "tsv", fields = "accession", size=1))
    if(status_code(res) == 200) {
      content_txt <- content(res, "text", encoding = "UTF-8")
      # TSV 返回格式通常有表头，第二行是 ID
      lines <- strsplit(content_txt, "\n")[[1]]
      if(length(lines) >= 2) {
        return(lines[2]) # 返回 Accession (例如 Q2M2I8)
      }
    }
    return(NULL)
  }, error = function(e) NULL)
}

cat(">>> 正在解析 ID 并尝试将 基因名 映射为 Uniprot ID (这需要一点时间)...\n")
site_list <- list()
pb <- txtProgressBar(min=0, max=length(unique_ids), style=3)

# 建立缓存避免重复查询同一个蛋白
gene_map_cache <- list()

for(i in seq_along(unique_ids)) {
  id_str <- unique_ids[i]
  # 解析原始 ID (假设格式 Gene_Pos_Mod)
  info <- tryCatch(parse_id_info(id_str), error=function(e) NULL)
  
  if(!is.null(info)) {
    gene_name <- info$acc # 这里目前解析出来的是 Gene Name (如 AAK1)
    pos <- as.numeric(info$pos)
    
    # 查找或使用缓存的 Uniprot ID
    real_acc <- NULL
    if(gene_name %in% names(gene_map_cache)) {
      real_acc <- gene_map_cache[[gene_name]]
    } else {
      # 这是一个 Gene Name，需要去查 ID
      mapped_id <- get_accession_from_gene(gene_name)
      if(!is.null(mapped_id)) {
        real_acc <- mapped_id
        gene_map_cache[[gene_name]] <- real_acc
      }
    }
    
    if(!is.null(real_acc)) {
      site_list[[i]] <- data.frame(
        Original_ID = id_str,
        Gene = gene_name,
        Uniprot = real_acc, # 这里存真正的 ID
        Position = pos,
        stringsAsFactors = F
      )
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

site_info <- bind_rows(site_list)
cat(sprintf("\n>>> 成功映射 %d 个位点 (共涉及 %d 个独立 Uniprot ID)。\n", nrow(site_info), length(unique(site_info$Uniprot))))

# --- 3. 批量下载 FASTA (使用真正的 Uniprot ID) ---
unique_proteins <- unique(site_info$Uniprot)
fasta_list <- list()
cat(sprintf(">>> 开始下载 %d 个蛋白序列...\n", length(unique_proteins)))
pb <- txtProgressBar(min=0, max=length(unique_proteins), style=3)

for(i in seq_along(unique_proteins)) {
  acc <- unique_proteins[i]
  url <- paste0("https://rest.uniprot.org/uniprotkb/", acc, ".fasta")
  
  tryCatch({
    res <- GET(url, config(connecttimeout = 10))
    if(status_code(res) == 200) {
      content_txt <- content(res, "text", encoding = "UTF-8")
      lines <- strsplit(content_txt, "\n")[[1]]
      # 过滤掉以 > 开头的行
      seq_lines <- lines[!grepl("^>", lines)]
      seq_str <- paste(seq_lines, collapse = "")
      
      if(nchar(seq_str) > 0) {
        fasta_list[[acc]] <- seq_str
      }
    }
  }, error = function(e) {})
  
  setTxtProgressBar(pb, i)
}
close(pb)

# --- 4. 生成 DeepMVP 文件 (修复空对象错误) ---

# [文件 A] FASTA
fasta_outfile <- "DeepMVP_Upload_Sequences.fasta"
if(length(fasta_list) > 0) {
  con <- file(fasta_outfile, "w")
  for(acc in names(fasta_list)) {
    writeLines(paste0(">", acc), con)
    writeLines(fasta_list[[acc]], con)
  }
  close(con)
  cat(sprintf("\n>>> FASTA 文件已生成: %s\n", fasta_outfile))
} else {
  stop("无法生成 FASTA 文件，下载列表为空！")
}

# [文件 B] Mutation List
mutation_outfile <- "DeepMVP_Batch_Mutation_Input.txt"
mut_lines <- character() # 初始化为字符向量

for(i in 1:nrow(site_info)) {
  acc <- site_info$Uniprot[i]
  pos <- site_info$Position[i]
  seq <- fasta_list[[acc]]
  
  if(!is.null(seq) && pos <= nchar(seq)) {
    ref_aa <- substr(seq, pos, pos)
    var_aa <- case_when(
      ref_aa == "S" ~ "A",
      ref_aa == "T" ~ "A",
      ref_aa == "Y" ~ "F",
      TRUE ~ "A"
    )
    # 格式: ID Ref Pos Var
    line <- sprintf("%s %s %d %s", acc, ref_aa, pos, var_aa)
    mut_lines <- c(mut_lines, line)
  }
}

if(length(mut_lines) > 0) {
  writeLines(mut_lines, mutation_outfile)
  cat(sprintf(">>> 突变列表已生成: %s (共 %d 条)\n", mutation_outfile, length(mut_lines)))
} else {
  cat("[Warning] 未生成突变数据。可能是位点位置与 Uniprot 最新序列不匹配。\n")
}