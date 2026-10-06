#!/usr/bin/env Rscript
## =============================================================================
## 09_gate1b 任务 2（Q2c）：作者 DAPAR 流程分支 —— 主流程（R）
##
## 运行（advisor 建的独立 R 环境，不碰全局 Python）：
##   <scratchpad>/envs/rdapar/bin/Rscript /home/user/phd-script/analysis/09_gate1b/02_author_branch.R
## 然后：
##   python3 /home/user/phd-script/analysis/09_gate1b/02_author_compare.py
##   （写 02_author_compare.tsv，并把 md 的 ⑥ 节填进 02_author_branch.md）
##
## 只读 data/、code/；输出全部写 analysis/09_gate1b/。参数全部按作者代码固定，
## 唯一的非作者值是 KNN 的 rng.seed = 1（作者为 sample(seq_len(1000), 1)，随机）。
## 作者脚本：code/Protein contour/Phospho/SpatialProteoDynamics.github.io-v1.0/
##   SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac/DataProcessing/
##   DAPAR_script_OSM.R（下称 D:）、translocation_plots.R（下称 T:）
## DAPAR 1.38 源码 deparse：<scratchpad>/dapar_1.38_src.txt（下称 src:）
## =============================================================================

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(limma)
  library(impute)
  library(metap)
  library(Biobase)
})

## ---------------------------------------------------------------- 路径
REPO  <- "/home/user/phd-script"
OUT   <- file.path(REPO, "analysis", "09_gate1b")
IN_PG <- file.path(REPO, "data", "pilot", "proteome_pg_matrix_120.tsv")
IN_PG_REL <- "data/pilot/proteome_pg_matrix_120.tsv"
AUTH_DIR <- file.path(REPO, "code", "Protein contour", "Phospho",
                      "SpatialProteoDynamics.github.io-v1.0",
                      "SpatialProteoDynamics-SpatialProteoDynamics.github.io-a6b8aac",
                      "DataProcessing")
AUTH_D <- file.path(AUTH_DIR, "DAPAR_script_OSM.R")
AUTH_T <- file.path(AUTH_DIR, "translocation_plots.R")
SCRATCH <- Sys.getenv("GATE1B_SCRATCH",
  "/tmp/claude-0/-home-user-phd-script/5ff66668-3d80-5c06-b5f4-becba97c6cd3/scratchpad")
DAPAR_SRC <- file.path(SCRATCH, "dapar_1.38_src.txt")

## ---------------------------------------------------------------- 参数（作者值）
TH       <- 3       # D:42  mvFilter(prot, type="atLeastOneCond", th=3)
SPAN     <- 0.7     # src:4 LOESS(..., span = 0.7)（D:67 未改默认）
KNN_K    <- 15      # D:68  wrapper.impute.KNN(prot_norm, 15)
ROWMAX   <- 0.99    # src:64
COLMAX   <- 0.99    # src:64
MAXP     <- 1500    # src:64
RNG_SEED <- 1       # src:65 作者为 sample(seq_len(1000), 1)；本轮固定为 1
QVAL     <- 0.025   # src:108 getQuantile4Imp 默认（D:72 未改）
FACTOR   <- 1       # src:108
TPS   <- c("2min", "8min", "20min", "90min")
CONDS <- c(TPS, "CTRL")          # 作者列布局：4 个处理条件块在前，CTRL 在 17:20（D:45、D:77）
FRS   <- paste0("FR", 1:6)
REPS  <- paste0("Rep", 1:4)
GENES4 <- c("GRB2", "SHC1", "CBL", "EGFR")
TP_ORDER_OUT <- c("CTRL", TPS)

## ---------------------------------------------------------------- 输出文件
F_BRANCH <- file.path(OUT, "02_author_branch.tsv")
F_CELLS  <- file.path(OUT, "02_author_cells.tsv")
F_MERGED <- file.path(OUT, "02_author_merged_log2.tsv")
F_DETQ   <- file.path(OUT, "02_author_detquant_values.tsv")
F_MD     <- file.path(OUT, "02_author_branch.md")
f_limma  <- function(j, tp) file.path(OUT, sprintf("02_author_limma_FR%d_%svsCTRL.tsv", j, tp))
R_OUTPUTS <- c(F_BRANCH, F_CELLS, F_MERGED, F_DETQ,
               unlist(lapply(1:6, function(j) vapply(TPS, function(tp) f_limma(j, tp), ""))))

## 复跑一致性：记录上一次运行的 md5（md 与 Python 产物除外）
prev_md5 <- if (all(file.exists(R_OUTPUTS))) tools::md5sum(R_OUTPUTS) else NULL

## ---------------------------------------------------------------- 工具函数
fmtn <- function(x) ifelse(is.na(x), "", sprintf("%.15g", x))
f4   <- function(x) ifelse(is.na(x), "", sprintf("%.4f", x))
g3   <- function(x) ifelse(is.na(x), "", sprintf("%.3g", x))
yn   <- function(b) ifelse(is.na(b), "", ifelse(b, "yes", "no"))
write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "", fileEncoding = "UTF-8")
}
md_table <- function(df) {
  df[] <- lapply(df, as.character)
  esc <- function(v) gsub("|", "\\|", as.character(v), fixed = TRUE)
  hdr <- paste0("| ", paste(esc(colnames(df)), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |"))
  c(hdr, sep, rows)
}
SELF <- data.frame(check = character(0), result = character(0), value = character(0))
add_check <- function(name, ok, value = "") {
  SELF[nrow(SELF) + 1, ] <<- list(name, if (isTRUE(ok)) "PASS" else "FAIL", as.character(value))
  if (!isTRUE(ok)) message("SELF-CHECK FAIL: ", name, " ", value)
}

## ---------------------------------------------------------------- 引用行号核对
## （验收 1：处理顺序表每步带作者行号与 DAPAR 源码行号 —— 这里逐条核对引用的行确实含该语句）
cite_check <- function(path, items) {
  L <- readLines(path, warn = FALSE)
  bad <- character(0)
  for (it in items) {
    ln <- it[[1]]; pat <- it[[2]]
    if (ln > length(L) || !grepl(pat, L[ln], fixed = TRUE)) bad <- c(bad, sprintf(":%d", ln))
  }
  list(n = length(items), bad = bad)
}
cites_D <- list(list(10, 'PG.Genes!=""'), list(11, 'gsub("^;", ""'), list(12, 'gsub(";.*", ""'),
  list(13, "!duplicated(exprsData$PG.Genes)"), list(21, 'fractions<-c("FR1"'),
  list(24, "contains(fractions[i])"), list(26, "grepl(fractions[i], Sample.name)"),
  list(40, "for(j in 1:6){"), list(41, "logData=T"),
  list(42, 'mvFilter(prot, type="atLeastOneCond", th=3)'),
  list(44, "length(x[!is.na(x)])"), list(48, "mv[i,5]>=3 || mv[i,1]>=3"),
  list(51, "mv[i,5]>=3 || mv[i,2]>=3"), list(54, "mv[i,5]>=3 || mv[i,3]>=3"),
  list(57, "mv[i,5]>=3 || mv[i,4]>=3"),
  list(67, 'wrapper.normalizeD(prot_mv, method="LOESS", type="overall")'),
  list(68, "wrapper.impute.KNN(prot_norm, 15)"), list(69, "findMECBlock(prot_norm)"),
  list(70, "reIntroduceMEC(prot_imp_KNN, MEC)"), list(72, "getQuantile4Imp(qData)$shiftedImpVal"),
  list(73, "impute.detQuant(qData, values)"), list(77, "prot_imp_2[,(4*k-3):(k*4)],prot_imp_2[,17:20]"),
  list(78, 'mv_keep[,k]=="KEEP"'),
  list(80, 'limmaCompleteTest(input_limma, conditions, comp.type="OnevsOne")'),
  list(81, 'method = "BH"'), list(89, "tapply(x,m, mean)"),
  list(93, "merge(merge_table, prot_imp_2, by=0, all=T)"),
  list(105, "rowSums(is.na(merge_table))<120"), list(107, "write.table(merge_table_na"))
cites_T <- list(list(7, "read.table(\"merged_table"), list(16, 'contains("2min")'),
  list(17, 'contains("CTRL")'), list(19, "2^prot_2min"), list(20, "2^prot_ctrl"),
  list(23, 'rep("FR1",4)'), list(24, "tapply(x,FR,function(x) {mean(x)})"),
  list(30, "x/sum(x)"), list(31, "x/sum(x)"),
  list(33, "abs(prot_2min_scaled-prot_Ctrl_scaled)"), list(34, "mean(x)"),
  list(36, "match(max(x),x)"), list(44, "sort(x,partial=len-N+1)[len-N+1]"),
  list(47, "match(maxN(x),x)"), list(56, "full_table$Gene.names"),
  list(63, "pval_FR1[pval_FR1$X==row.names(temp_table)[i],3]"), list(65, "pval[i]=1"),
  list(152, "sumlog(c(pval[i],pval_2[i]))$p"), list(159, "limegreen"), list(163, "dodgerblue"),
  list(167, "coral"), list(170, 'color[i]="gray"'),
  list(173, 'p.adjust(temp_table$pval_combi, method="BH")'),
  list(181, "(temp_table$pval_combi_FDR<0.05) & (temp_table$MS_2min_mean>0.1)"))
cites_S <- list(list(11, 'limma::normalizeCyclicLoess(x = qData, method = "fast"'),
  list(12, "span = span"), list(43, "LOESS(qData, ...)"),
  list(59, "conditions <- unique(Biobase::pData(obj)$Condition)"),
  list(63, "impute::impute.knn(Biobase::exprs(obj)[, ind]"),
  list(64, "k = K, rowmax = 0.99, colmax = 0.99, maxp = 1500"),
  list(65, "rng.seed = sample(seq_len(1000), 1)"), list(68, "== 0] <- NA"),
  list(83, "lNA <- which(apply(is.na(Biobase::exprs(obj)[, ind]"),
  list(102, "as.vector(replicates)] <- NA"), list(108, "qval = 0.025, factor = 1"),
  list(111, "apply(qdata, 2, stats::quantile, qval, na.rm = TRUE)"), list(112, "r1 * factor"),
  list(169, "limma::eBayes(limma::contrasts.fit(limma::lmFit(qData"),
  list(217, "factor(sTab$Condition, ordered = TRUE)"), list(234, "design[seq, j] <- rep(1, length(seq))"),
  list(281, 'paste("(", label.agg[i], ")/"'))
cc <- list(D = cite_check(AUTH_D, cites_D), T = cite_check(AUTH_T, cites_T),
           S = if (file.exists(DAPAR_SRC)) cite_check(DAPAR_SRC, cites_S) else list(n = 0, bad = "缺 dapar_1.38_src.txt"))
add_check("引用行号与作者脚本 / DAPAR 源码逐行一致",
          length(cc$D$bad) + length(cc$T$bad) + length(cc$S$bad) == 0,
          sprintf("DAPAR_script_OSM.R %d 处、translocation_plots.R %d 处、dapar_1.38_src.txt %d 处；不一致 %s",
                  cc$D$n, cc$T$n, cc$S$n,
                  if (length(c(cc$D$bad, cc$T$bad, cc$S$bad))) paste(c(cc$D$bad, cc$T$bad, cc$S$bad), collapse = ",") else "0"))

## ---------------------------------------------------------------- DAPAR 1.38（lazyLoad，仅作交叉核对与 make.design.1 / make.contrast）
dapar <- new.env()
invisible(lazyLoad(file.path(.libPaths()[1], "DAPAR", "R", "DAPAR"), envir = dapar))
for (nm in ls(dapar)) if (is.function(dapar[[nm]])) environment(dapar[[nm]]) <- dapar
## （lazyLoad 取出的闭包环境为 globalenv，找不到 pkgs.require 等内部函数，故把闭包环境改为该 env）

## ---------------------------------------------------------------- 读表
raw <- read.delim(IN_PG, sep = "\t", colClasses = "character", na.strings = character(0),
                  quote = "", check.names = FALSE, comment.char = "")
stopifnot(nrow(raw) == 10000, ncol(raw) == 122,
          identical(colnames(raw)[1:2], c("Protein.Group", "Genes")))

## run 名解析：07 parse_design 的三个正则
parse_design <- function(run) {
  m1 <- regmatches(run, regexec("_(2min|8min|20min|90min|CTRL)_", run, ignore.case = TRUE))[[1]]
  m2 <- regmatches(run, regexec("_(FR[0-9])_", run))[[1]]
  m3 <- regmatches(run, regexec("_(Rep[0-9])", run))[[1]]
  if (length(m1) < 2 || length(m2) < 2 || length(m3) < 2) return(c(NA, NA, NA))
  c(m1[2], m2[2], m3[2])
}
runs <- colnames(raw)[-(1:2)]
dz <- t(vapply(runs, parse_design, character(3)))
design <- data.frame(run = runs, timepoint = dz[, 1], fraction = dz[, 2], rep = dz[, 3],
                     row.names = NULL)
stopifnot(!anyNA(design), all(design$timepoint %in% CONDS), all(design$fraction %in% FRS),
          all(design$rep %in% REPS), nrow(unique(design[, 2:4])) == 120)
run_of <- function(tp, fr, rp) design$run[design$timepoint == tp & design$fraction == fr & design$rep == rp]

## ---------------------------------------------------------------- S0 去空基因名、取第一个基因名、去重（D:10-13）
n_raw <- nrow(raw)
s0 <- raw[raw$Genes != "", , drop = FALSE]                    # D:10
n_drop_empty <- n_raw - nrow(s0)
s0$Genes <- gsub("^;", "", s0$Genes)                          # D:11
s0$Genes <- gsub(";.*", "", s0$Genes)                         # D:12
n_multi_gene <- sum(grepl(";", raw$Genes[raw$Genes != ""]))
dup <- duplicated(s0$Genes)
n_drop_dup <- sum(dup)
s0 <- s0[!dup, , drop = FALSE]                                # D:13
genes <- s0$Genes
stopifnot(!anyDuplicated(genes), all(GENES4 %in% genes))
n_s0 <- nrow(s0)
pg_of4 <- setNames(s0$Protein.Group[match(GENES4, genes)], GENES4)

## 线性值矩阵（空串 = 缺失）
Xlin <- vapply(runs, function(r) {
  v <- s0[[r]]; out <- rep(NA_real_, length(v)); ok <- v != ""
  out[ok] <- suppressWarnings(as.numeric(v[ok])); out
}, numeric(n_s0))
rownames(Xlin) <- genes
n_nonempty <- sum(vapply(runs, function(r) sum(s0[[r]] != ""), 0))
n_parse_fail <- n_nonempty - sum(!is.na(Xlin))
stopifnot(n_parse_fail == 0)
raw_nonempty4 <- vapply(GENES4, function(g) sum(!is.na(Xlin[g, ])), 0)

## ---------------------------------------------------------------- S2 log2（D:41 createMSnset(..., logData=T)）
n_nonpos <- sum(Xlin <= 0, na.rm = TRUE)
L2 <- Xlin
L2[!is.na(L2) & L2 <= 0] <- NA
L2 <- log2(L2)

## ---------------------------------------------------------------- impute.knn 的逐行复刻（只为计数）
## 复刻 impute 1.80.0 的 impute.knn → knnimp → knnimp.internal / knnimp.split，
## 记录：缺失比例 > rowmax 的行（列均值填补）、KNN 后仍为 NA 的行（列均值填补）、
## 被二分聚类切成 ≤ k 行的小簇里的缺失格（库代码对其调 meanimp(x[index, ])，但此时缺失格已被置 0，
## meanimp 不识别，结果留 0 —— 由 DAPAR src:68 置回 NA）。复刻输出必须与 impute::impute.knn 完全相同。
TR <- new.env()
knnimp_internal_tr <- function(x, k, imiss, irmiss, p, n, maxp) {
  if (p <= maxp) impute:::knnimp.internal(x, k, imiss, irmiss, p, n, maxp = maxp)
  else knnimp_split_tr(x, k, imiss, irmiss, p, n, maxp = maxp)
}
knnimp_split_tr <- function(x, k, imiss, irmiss, p, n, maxp) {
  junk <- impute:::twomeans.miss(x)
  size <- junk$size
  clus <- junk$cluster
  TR$n_split <- TR$n_split + 1L
  for (i in seq(size)) {
    p <- as.integer(size[i])
    index <- clus == i
    if (p <= k) {
      im <- imiss[index, , drop = FALSE]
      TR$small_rows  <- TR$small_rows + sum(rowSums(im) > 0)
      TR$small_cells <- TR$small_cells + sum(im)
    }
    x[index, ] <- if (p <= k) impute:::meanimp(x[index, ])
                  else knnimp_internal_tr(x[index, ], k, imiss[index, ], irmiss[index], p, n, maxp)
  }
  x
}
knnimp_tr <- function(x, k = 10, maxmiss = 0.5, maxp = 1500) {
  pn <- dim(x); dn <- dimnames(x)
  p <- as.integer(pn[1]); n <- as.integer(pn[2])
  imiss <- is.na(x); x[imiss] <- 0
  irmiss <- drop(imiss %*% rep(1, n))
  imax <- trunc(maxmiss * n); imax <- irmiss > imax; simax <- sum(imax)
  TR$simax <- simax
  if (simax > 0) {
    irmiss <- irmiss[!imax]; imiss.omit <- imiss[imax, , drop = FALSE]; imiss <- imiss[!imax, ]
    xomit <- x[imax, , drop = FALSE]; x <- x[!imax, ]; p <- as.integer(p - simax)
  }
  storage.mode(imiss) <- "integer"; storage.mode(irmiss) <- "integer"; storage.mode(x) <- "double"
  ximp <- if (p <= maxp) impute:::knnimp.internal(x, k, imiss, irmiss, p, n, maxp = maxp)
          else knnimp_split_tr(x, k, imiss, irmiss, p, n, maxp = maxp)
  imiss.new <- is.na(ximp); newmiss <- any(imiss.new)
  TR$newmiss_rows <- sum(rowSums(imiss.new) > 0)
  TR$newmiss_cells <- sum(imiss.new)
  nm_full <- matrix(FALSE, pn[1], pn[2])
  if (simax > 0) nm_full[!imax, ] <- imiss.new else nm_full[] <- imiss.new
  TR$newmiss_mask <- nm_full
  if ((simax > 0) | newmiss) {
    xbar <- impute:::mean.miss(x, imiss = imiss)
    if (newmiss) ximp <- impute:::meanimp(ximp, imiss.new, xbar)
    if (simax > 0) {
      xomit <- impute:::meanimp(xomit, imiss.omit, xbar)
      xout <- array(0, dim = pn); xout[!imax, ] <- ximp; xout[imax, ] <- xomit; ximp <- xout
    }
  }
  dimnames(ximp) <- dn
  ximp
}
impute_knn_tr <- function(data, k, rowmax, colmax, maxp, rng.seed) {
  for (v in c("simax", "newmiss_rows", "newmiss_cells", "small_rows", "small_cells", "n_split")) TR[[v]] <- 0L
  set.seed(rng.seed)
  x <- data; p <- nrow(x)
  col.nas <- drop(rep(1, p) %*% is.na(x))
  if (any(col.nas > colmax * p)) stop("a column has too many missing values")
  out <- knnimp_tr(x, k, maxmiss = rowmax, maxp = maxp)
  list(data = out, simax = TR$simax, newmiss_rows = TR$newmiss_rows, newmiss_cells = TR$newmiss_cells,
       newmiss_mask = TR$newmiss_mask,
       small_rows = TR$small_rows, small_cells = TR$small_cells, n_split = TR$n_split)
}

mk_es <- function(m, cond) {
  pd <- data.frame(Sample.name = colnames(m), Condition = cond, Bio.Rep = seq_along(cond),
                   row.names = colnames(m))
  Biobase::ExpressionSet(assayData = m, phenoData = Biobase::AnnotatedDataFrame(pd))
}

## ---------------------------------------------------------------- 每个 fraction：S1–S8
FRSTAT <- list(); LIMMA <- list(); FINAL <- list(); SRC <- list(); DETQ <- list(); MVK <- list()
COND_STAT <- list(); XCHK <- list(); COLMEAN <- list()
for (j in 1:6) {
  fr <- FRS[j]
  ## S1（D:21-28、D:40）：该 fraction 的 20 列，按 [2min, 8min, 20min, 90min, CTRL] × Rep1..4 排
  cols <- unlist(lapply(CONDS, function(cn) vapply(REPS, function(rp) run_of(cn, fr, rp), "")))
  stopifnot(length(cols) == 20, !anyDuplicated(cols))
  cond <- rep(CONDS, each = 4)
  Lj <- L2[, cols, drop = FALSE]

  ## S3 mvFilter(type="atLeastOneCond", th=3)（D:42；旧版逻辑：至少一个条件非缺失数 ≥ th）
  cnt <- sapply(CONDS, function(cn) rowSums(!is.na(Lj[, cond == cn, drop = FALSE])))
  keep <- apply(cnt >= TH, 1, any)
  Xf <- Lj[keep, , drop = FALSE]
  nk <- nrow(Xf)

  ## S4 对比保留标记（D:44-58）：CTRL ≥ 3 或 条件 k ≥ 3
  mv <- cnt[keep, , drop = FALSE]
  KEEP <- sapply(TPS, function(tp) mv[, "CTRL"] >= 3 | mv[, tp] >= 3)
  MVK[[fr]] <- list(mv = mv, KEEP = KEEP)

  ## S5 LOESS（D:67 → src:42-43 → src:10-13）
  N <- limma::normalizeCyclicLoess(x = Xf, method = "fast", span = SPAN)
  x_loess <- identical(N, dapar$LOESS(Xf, conds = cond, type = "overall"))
  stopifnot(identical(is.na(N), is.na(Xf)))

  ## S6 KNN（D:68 → src:59-67），逐条件 4 列，rng.seed = 1
  Kimp <- N
  colmean_mask <- matrix(FALSE, nk, 20, dimnames = dimnames(Xf))   # impute.knn 对 KNN 后仍为 NA 的格退回列均值
  cstat <- data.frame()
  knn_identical <- TRUE
  for (cn in unique(cond)) {
    ind <- which(cond == cn)
    wmsg <- character(0)
    co <- capture.output(
      res <- withCallingHandlers(
        impute::impute.knn(N[, ind], k = KNN_K, rowmax = ROWMAX, colmax = COLMAX,
                           maxp = MAXP, rng.seed = RNG_SEED),
        warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") }))
    tr <- impute_knn_tr(N[, ind], k = KNN_K, rowmax = ROWMAX, colmax = COLMAX, maxp = MAXP,
                        rng.seed = RNG_SEED)
    knn_identical <- knn_identical && identical(res$data, tr$data)
    warn_rows <- if (length(wmsg)) sum(as.integer(sub("^([0-9]+) rows with more than.*", "\\1",
                                                    wmsg[grepl("rows with more than", wmsg)]))) else 0L
    Kimp[, ind] <- res$data
    colmean_mask[, ind] <- tr$newmiss_mask
    cstat <- rbind(cstat, data.frame(fraction = fr, condition = cn,
      missing_cells = sum(is.na(N[, ind])),
      knn_meanimp_rows_rowmax = tr$simax, knn_warning_rows = warn_rows,
      knn_meanimp_rows_newmiss = tr$newmiss_rows, knn_meanimp_cells_newmiss = tr$newmiss_cells,
      knn_smallcluster_rows = tr$small_rows, knn_smallcluster_cells = tr$small_cells,
      knn_cluster_lines = sum(grepl("^Cluster size", co)),
      MEC_blocks = sum(rowSums(is.na(N[, ind])) == length(ind))))
  }
  zero_mask <- !is.na(Kimp) & Kimp == 0
  n_zero <- sum(zero_mask)
  Kimp[zero_mask] <- NA                                        # src:68

  ## S7a MEC（D:69-70 → src:76-93、src:96-105）
  ucond <- unique(cond)
  MEC <- data.frame()
  for (ci in seq_along(ucond)) {
    ind <- which(cond == ucond[ci])
    lNA <- which(apply(is.na(N[, ind]), 1, sum) == length(ind))
    if (length(lNA) > 0) {
      tmp <- data.frame(ci, lNA); names(tmp) <- c("Condition", "Line"); MEC <- rbind(MEC, tmp)
    }
  }
  M <- Kimp
  for (i in seq_len(nrow(MEC))) M[MEC[i, "Line"], which(cond == ucond[MEC[i, "Condition"]])] <- NA
  mec_mask <- matrix(FALSE, nk, 20)
  for (i in seq_len(nrow(MEC))) mec_mask[MEC[i, "Line"], which(cond == ucond[MEC[i, "Condition"]])] <- TRUE
  ## DAPAR 1.38 原函数交叉核对（ExpressionSet 包装）
  MEC_d <- dapar$findMECBlock(mk_es(N, cond))
  x_mec <- identical(unname(as.matrix(MEC_d)), unname(as.matrix(MEC)))
  M_d <- Biobase::exprs(dapar$reIntroduceMEC(mk_es(Kimp, cond), MEC_d))
  x_reint <- identical(unname(M_d), unname(M))

  ## S7b 低值填补：getQuantile4Imp(qval=0.025, factor=1)（D:72 → src:108-114）+ 旧版 impute.detQuant（D:73）
  qv <- apply(M, 2, stats::quantile, QVAL, na.rm = TRUE)       # src:111
  shifted <- qv * FACTOR                                        # src:112
  x_q <- isTRUE(all.equal(unname(dapar$getQuantile4Imp(M)$shiftedImpVal), unname(shifted), tolerance = 0))
  Fm <- M
  for (i in seq_len(ncol(Fm))) {                                # 旧版 impute.detQuant(qData, values)
    col <- Fm[, i]; col[which(is.na(col))] <- shifted[i]; Fm[, i] <- col
  }
  stopifnot(!anyNA(Fm))

  ## 每格来源
  src <- matrix("observed", nk, 20, dimnames = dimnames(Xf))
  src[is.na(Xf)] <- "KNN"
  src[is.na(M)] <- "detQuant"
  stopifnot(all(is.na(Xf[is.na(M)])))                         # detQuant 只落在原缺失格
  stopifnot(all(src[colmean_mask] == "KNN"))                   # 列均值退回格不在 MEC 块内，来源记 KNN
  COLMEAN[[fr]] <- colmean_mask
  obs_equal <- max(abs(Fm[!is.na(Xf)] - N[!is.na(Xf)]))         # observed 格 = LOESS 后观测值
  DETQ[[fr]] <- data.frame(run = cols, fraction = fr, timepoint = cond, rep = rep(REPS, 5),
                           q025_log2 = unname(qv), shiftedImpVal_log2 = unname(shifted),
                           n_values_used = colSums(!is.na(M)), n_detQuant_cells = colSums(src == "detQuant"))

  n_obs <- sum(src == "observed"); n_knn <- sum(src == "KNN"); n_det <- sum(src == "detQuant")
  FRSTAT[[fr]] <- data.frame(fraction = fr, rows_in = nrow(Lj), S3_kept = nk,
    observed = n_obs, KNN = n_knn, MEC_blocks = nrow(MEC), detQuant = n_det,
    detQuant_MEC = sum(mec_mask), detQuant_zero = sum(zero_mask & !mec_mask),
    knn_meanimp_rows = sum(cstat$knn_meanimp_rows_rowmax) + sum(cstat$knn_meanimp_rows_newmiss),
    knn_meanimp_rows_rowmax = sum(cstat$knn_meanimp_rows_rowmax),
    knn_meanimp_rows_newmiss = sum(cstat$knn_meanimp_rows_newmiss),
    knn_meanimp_cells_newmiss = sum(colmean_mask),
    knn_warning_rows = sum(cstat$knn_warning_rows),
    knn_zero_cells = n_zero, knn_smallcluster_cells = sum(cstat$knn_smallcluster_cells),
    sum_ok = (n_obs + n_knn + n_det) == nk * 20,
    KEEP_2min = sum(KEEP[, "2min"]), KEEP_8min = sum(KEEP[, "8min"]),
    KEEP_20min = sum(KEEP[, "20min"]), KEEP_90min = sum(KEEP[, "90min"]))
  COND_STAT[[fr]] <- cstat

  ## S8 limma（D:76-82 → src:146-194）：非配对，条件指示矩阵（src:215-239），BH（D:81）
  sTab_full <- data.frame(Sample.name = cols, Condition = cond, Bio.Rep = seq_along(cols))
  lim_maxdiff <- 0
  for (k in 1:4) {
    tp <- TPS[k]
    idx <- c(which(cond == tp), which(cond == "CTRL"))         # D:77 (4k-3):(4k) 与 17:20
    inp <- Fm[, idx, drop = FALSE][KEEP[, tp], , drop = FALSE] # D:78
    sTab <- sTab_full[idx, ]                                   # D:79
    dm <- dapar$make.design.1(sTab)                            # src:215-239
    conds_f <- factor(sTab$Condition, levels = unique(sTab$Condition))   # src:160
    contra <- dapar$make.contrast(dm, condition = conds_f, contrast = 1, design.level = 1)  # src:242-301
    cmtx <- limma::makeContrasts(contrasts = contra, levels = make.names(colnames(dm)))      # src:167-168
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(inp, dm), cmtx))                # src:169-170
    logFC <- unname(fit$coefficients[, 1]); pv <- unname(fit$p.value[, 1])
    bh <- p.adjust(pv, method = "BH", n = length(pv))          # D:81
    ref <- dapar$limmaCompleteTest(inp, sTab, comp.type = "OnevsOne")
    lim_maxdiff <- max(lim_maxdiff, abs(ref$logFC[, 1] - logFC), abs(ref$P_Value[, 1] - pv))
    stopifnot(identical(rownames(ref$logFC), rownames(inp)),
              grepl(paste0("^", tp, "_vs_CTRL"), colnames(ref$logFC)[1]))
    LIMMA[[paste(fr, tp)]] <- data.frame(gene = rownames(inp), logFC = logFC, P_Value = pv, BH = bh)
  }
  XCHK[[fr]] <- data.frame(fraction = fr, loess_identical = x_loess, knn_trace_identical = knn_identical,
                           findMEC_identical = x_mec, reIntroduceMEC_identical = x_reint,
                           getQuantile4Imp_identical = x_q, limma_vs_limmaCompleteTest_maxdiff = lim_maxdiff,
                           observed_vs_loess_maxdiff = obs_equal)
  FINAL[[fr]] <- Fm; SRC[[fr]] <- src
  cat(sprintf("%s kept=%d obs=%d KNN=%d MEC=%d detQ=%d zero=%d\n", fr, nk, n_obs, n_knn, nrow(MEC), n_det, n_zero))
}
FRSTAT <- do.call(rbind, FRSTAT); COND_STAT <- do.call(rbind, COND_STAT); XCHK <- do.call(rbind, XCHK)
DETQ <- do.call(rbind, DETQ)

## ---------------------------------------------------------------- 写 limma 表（24 个）
for (j in 1:6) for (tp in TPS) {
  t0 <- LIMMA[[paste(FRS[j], tp)]]
  write_tsv(data.frame(gene = t0$gene, logFC = fmtn(t0$logFC), P_Value = fmtn(t0$P_Value),
                       BH = fmtn(t0$BH)), f_limma(j, tp))
}

## ---------------------------------------------------------------- S9 合并（D:89-107）
mcols <- unlist(lapply(FRS, function(fr) unlist(lapply(CONDS, function(cn)
  vapply(REPS, function(rp) run_of(cn, fr, rp), "")))))
Fall <- matrix(NA_real_, n_s0, 120, dimnames = list(genes, mcols))
Sall <- matrix("filtered_out", n_s0, 120, dimnames = list(genes, mcols))
for (fr in FRS) {
  Fm <- FINAL[[fr]]
  Fall[rownames(Fm), colnames(Fm)] <- Fm
  Sall[rownames(Fm), colnames(Fm)] <- SRC[[fr]]
}
merged_na <- Fall[rowSums(is.na(Fall)) < 120, , drop = FALSE]    # D:105
n_pass_tab <- table(factor(rowSums(sapply(FRS, function(fr) genes %in% rownames(FINAL[[fr]]))), levels = 0:6))
n_merged <- nrow(merged_na)
mdf <- data.frame(gene = rownames(merged_na), check.names = FALSE)
for (cname in mcols) mdf[[cname]] <- fmtn(merged_na[, cname])
write_tsv(mdf, F_MERGED)

## detQuant 值
write_tsv(data.frame(run = DETQ$run, fraction = DETQ$fraction, timepoint = DETQ$timepoint, rep = DETQ$rep,
                     q025_log2 = fmtn(DETQ$q025_log2), shiftedImpVal_log2 = fmtn(DETQ$shiftedImpVal_log2),
                     n_values_used = DETQ$n_values_used, n_detQuant_cells = DETQ$n_detQuant_cells), F_DETQ)

## ---------------------------------------------------------------- 迁移分析 T1–T10（每个 EGF 时间点一遍）
maxN <- function(x, N = 2) {                                     # T:38-45 原样
  len <- length(x)
  if (N > len) {
    warning("N greater than length(x).  Setting N=length(x)")
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}
## 作者 T:47 的 maxN 对含 NA 的行报错（sort.int: index 5 outside bounds），作者代码只能在
## 合并表无 NA 的行上运行；这里把 T1–T10 限于 6 个 fraction 都过 S3 的蛋白（120 列全非 NA）。
complete <- rowSums(is.na(merged_na)) == 0
TT <- merged_na[complete, , drop = FALSE]
n_complete <- nrow(TT)
maxN_err_demo <- inherits(try(maxN(rep(NA_real_, 6)), silent = TRUE), "try-error")
color_of <- function(m1, m2) {                                   # T:157-171
  if ((m1 == 1 || m1 == 2) && (m2 == 3 || m2 == 4)) "limegreen"
  else if ((m1 == 3 || m1 == 4) && (m2 == 1 || m2 == 2)) "limegreen"
  else if ((m1 == 1 || m1 == 2) && (m2 == 5 || m2 == 6)) "dodgerblue"
  else if ((m1 == 5 || m1 == 6) && (m2 == 1 || m2 == 2)) "dodgerblue"
  else if ((m1 == 3 || m1 == 4) && (m2 == 5 || m2 == 6)) "coral"
  else if ((m1 == 5 || m1 == 6) && (m2 == 3 || m2 == 4)) "coral"
  else "gray"
}
FRlab <- rep(FRS, each = 4)                                       # T:23
TRANS <- list(); PROT_SUM <- data.frame()
for (tp in TPS) {
  ctp <- unlist(lapply(FRS, function(fr) vapply(REPS, function(rp) run_of(tp, fr, rp), "")))
  cct <- unlist(lapply(FRS, function(fr) vapply(REPS, function(rp) run_of("CTRL", fr, rp), "")))
  A <- 2^TT[, ctp, drop = FALSE]; B <- 2^TT[, cct, drop = FALSE]  # T:16-20
  Am <- t(apply(A, 1, function(x) tapply(x, FRlab, function(x) { mean(x) })))   # T:24
  Bm <- t(apply(B, 1, function(x) tapply(x, FRlab, function(x) { mean(x) })))   # T:25
  As <- t(apply(Am, 1, function(x) { x / sum(x) }))                             # T:30
  Bs <- t(apply(Bm, 1, function(x) { x / sum(x) }))                             # T:31
  MS <- abs(As - Bs)                                                            # T:33
  MSmean <- apply(MS, 1, function(x) mean(x))                                   # T:34
  max1 <- apply(MS, 1, function(x) { match(max(x), x) })                        # T:36
  max2 <- apply(MS, 1, function(x) { match(maxN(x), x) })                       # T:47
  ## T:58-148：fraction max1 / max2 的 limma 表第 3 列（P_Value）；蛋白不在表 → 1
  ptab <- lapply(FRS, function(fr) LIMMA[[paste(fr, tp)]])
  getp <- function(k, g) { t0 <- ptab[[k]]; x <- t0$P_Value[t0$gene == g]; if (length(x) == 0) 1 else x }
  gn <- rownames(TT)
  p1 <- mapply(getp, max1, gn); p2 <- mapply(getp, max2, gn)
  pc <- vapply(seq_along(p1), function(i) metap::sumlog(c(p1[i], p2[i]))$p, numeric(1))   # T:150-153
  colr <- mapply(color_of, max1, max2)                                          # T:157-171
  fdr <- p.adjust(pc, method = "BH")                                            # T:173
  TRANS[[tp]] <- list(As = As, Bs = Bs, MSmean = MSmean, max1 = max1, max2 = max2,
                      p1 = p1, p2 = p2, pc = pc, colr = colr, fdr = fdr, n = length(pc))
  PROT_SUM <- rbind(PROT_SUM, data.frame(timepoint = tp, n_proteins_in_FDR = sum(!is.na(pc)),
    n_MS_gt_0.1 = sum(MSmean > 0.1), n_FDR_lt_0.05 = sum(fdr < 0.05),
    n_label_T181 = sum(fdr < 0.05 & MSmean > 0.1),
    n_p1_eq_1 = sum(p1 == 1), n_p2_eq_1 = sum(p2 == 1)))
}

## ---------------------------------------------------------------- 02_author_branch.tsv（16 行）
br <- data.frame()
for (g in GENES4) for (tp in TPS) {
  R <- TRANS[[tp]]
  if (g %in% rownames(TT)) {
    m1 <- R$max1[[g]]; m2 <- R$max2[[g]]
    s_tp <- R$As[g, m1]; s_c <- R$Bs[g, m1]
    br <- rbind(br, data.frame(gene = g, timepoint = tp, MS_author_pipeline = fmtn(R$MSmean[[g]]),
      MS_max1 = FRS[m1], MS_max2 = FRS[m2], share_max1_tp = fmtn(s_tp), share_max1_ctrl = fmtn(s_c),
      direction_max1 = if (s_tp > s_c) "up" else if (s_tp < s_c) "down" else "equal",
      p_max1 = fmtn(R$p1[[g]]), p_max2 = fmtn(R$p2[[g]]),
      pval_combi = fmtn(R$pc[match(g, rownames(TT))]), pval_combi_FDR = fmtn(R$fdr[match(g, rownames(TT))]),
      color_category = R$colr[[g]], MS_gt_0.1 = yn(R$MSmean[[g]] > 0.1),
      FDR_lt_0.05 = yn(R$fdr[match(g, rownames(TT))] < 0.05), n_proteins_in_FDR = R$n, check.names = FALSE))
  } else {
    br <- rbind(br, data.frame(gene = g, timepoint = tp, MS_author_pipeline = "", MS_max1 = "", MS_max2 = "",
      share_max1_tp = "", share_max1_ctrl = "", direction_max1 = "", p_max1 = "", p_max2 = "",
      pval_combi = "", pval_combi_FDR = "", color_category = "", MS_gt_0.1 = "", FDR_lt_0.05 = "",
      n_proteins_in_FDR = R$n, check.names = FALSE))
  }
}
stopifnot(nrow(br) == 16)
write_tsv(br, F_BRANCH)

## ---------------------------------------------------------------- 02_author_cells.tsv（480 行）
cells <- data.frame()
for (g in GENES4) for (tp in TP_ORDER_OUT) for (fr in FRS) for (rp in REPS) {
  r <- run_of(tp, fr, rp)
  cells <- rbind(cells, data.frame(gene = g, timepoint = tp, rep = rp, fraction = fr,
                                   value_log2_final = fmtn(Fall[g, r]), source = Sall[g, r]))
}
stopifnot(nrow(cells) == 480)
write_tsv(cells, F_CELLS)

## ---------------------------------------------------------------- 自检（R 部分）
add_check("每 fraction：observed + KNN + detQuant = S3 保留行数 × 20", all(FRSTAT$sum_ok),
          paste(sprintf("%s %d+%d+%d=%d×20", FRSTAT$fraction, FRSTAT$observed, FRSTAT$KNN, FRSTAT$detQuant,
                        FRSTAT$S3_kept), collapse = "；"))
add_check("impute.knn 退回列均值的行（缺失 > rowmax）= MEC 块数（每 fraction）",
          all(FRSTAT$knn_meanimp_rows_rowmax == FRSTAT$MEC_blocks) &&
            all(COND_STAT$knn_warning_rows == COND_STAT$knn_meanimp_rows_rowmax),
          "与 impute.knn 警告中的行数一致")
add_check("detQuant 格 = MEC 格 + KNN 输出为 0 的格（src:68）",
          all(FRSTAT$detQuant == FRSTAT$detQuant_MEC + FRSTAT$detQuant_zero), "")
add_check("KNN_colmean 格（impute.knn 内 KNN 后仍为 NA → 列均值）全部来源记 KNN、格数 = 复刻计数",
          all(FRSTAT$knn_meanimp_cells_newmiss == aggregate(COND_STAT$knn_meanimp_cells_newmiss,
                                                           list(COND_STAT$fraction), sum)$x[match(FRSTAT$fraction, sort(unique(COND_STAT$fraction)))]),
          paste(sprintf("%s %d", FRSTAT$fraction, FRSTAT$knn_meanimp_cells_newmiss), collapse = " / "))
add_check("impute.knn 逐行复刻（计数用）与 impute::impute.knn 输出完全相同", all(XCHK$knn_trace_identical), "30 次调用")
add_check("S5 = DAPAR 1.38 LOESS(type='overall')（identical）", all(XCHK$loess_identical), "6 fraction")
add_check("S7a = DAPAR 1.38 findMECBlock / reIntroduceMEC（identical）",
          all(XCHK$findMEC_identical) && all(XCHK$reIntroduceMEC_identical), "6 fraction")
add_check("S7b 分位数 = DAPAR 1.38 getQuantile4Imp（identical）", all(XCHK$getQuantile4Imp_identical), "120 列")
add_check("S8 直接 limma = DAPAR 1.38 limmaCompleteTest（logFC、P_Value）",
          max(XCHK$limma_vs_limmaCompleteTest_maxdiff) < 1e-12,
          sprintf("max|差| = %.1e", max(XCHK$limma_vs_limmaCompleteTest_maxdiff)))
add_check("observed 格的最终值 = LOESS 后观测值", max(XCHK$observed_vs_loess_maxdiff) == 0,
          sprintf("max|差| = %.1e", max(XCHK$observed_vs_loess_maxdiff)))
## 四蛋白 observed 格数 = 原始非空格数 − 落在没过 S3 的 fraction 里的非空格数
g_obs <- vapply(GENES4, function(g) sum(cells$gene == g & cells$source == "observed"), 0)
g_fo_raw <- vapply(GENES4, function(g) {
  sum(vapply(FRS, function(fr) if (g %in% rownames(FINAL[[fr]])) 0 else
    sum(!is.na(Xlin[g, design$run[design$fraction == fr]])), 0))
}, 0)
add_check("四蛋白原始非空格数 = GRB2 90 / SHC1 115 / CBL 59 / EGFR 80",
          identical(unname(raw_nonempty4), c(90, 115, 59, 80)),
          paste(sprintf("%s %d", GENES4, raw_nonempty4), collapse = " / "))
add_check("四蛋白 observed 格数 = 原始非空格数 − 没过 S3 的 fraction 中的非空格数",
          all(g_obs == raw_nonempty4 - g_fo_raw),
          paste(sprintf("%s %d = %d − %d", GENES4, g_obs, raw_nonempty4, g_fo_raw), collapse = "；"))
## cells 与合并矩阵一致（同一 Fall/Sall 写出；这里按写出的字符串回读核对）
mchk <- read.delim(F_MERGED, sep = "\t", colClasses = "character", na.strings = character(0),
                   quote = "", check.names = FALSE)
cmis <- 0
for (i in seq_len(nrow(cells))) {
  r <- run_of(cells$timepoint[i], cells$fraction[i], cells$rep[i])
  w <- match(cells$gene[i], mchk$gene)
  mv0 <- if (is.na(w)) "" else mchk[[r]][w]
  if (!identical(mv0, cells$value_log2_final[i])) cmis <- cmis + 1
  if ((cells$source[i] == "filtered_out") != (cells$value_log2_final[i] == "")) cmis <- cmis + 1
}
add_check("02_author_cells 的值与 02_author_merged_log2 一致（filtered_out ⇔ 空）", cmis == 0,
          sprintf("480 格，不一致 %d", cmis))
add_check("branch 16 行；n_proteins_in_FDR = 合并表中 120 列全非 NA 的行数",
          nrow(br) == 16 && all(br$n_proteins_in_FDR == n_complete), sprintf("%d", n_complete))
add_check("pval_combi_FDR 的 BH 在全部可算蛋白上做（非四蛋白）", all(PROT_SUM$n_proteins_in_FDR == n_complete),
          paste(sprintf("%s %d", PROT_SUM$timepoint, PROT_SUM$n_proteins_in_FDR), collapse = " / "))
add_check("作者 maxN 对含 NA 的行报错（确认 T1–T10 限于完整行的原因）", maxN_err_demo, "sort.int: index 5 outside bounds")
new_md5 <- tools::md5sum(R_OUTPUTS)
rerun_txt <- if (is.null(prev_md5)) "无上一次运行的产物（首次运行）" else {
  same <- sum(unname(prev_md5) == unname(new_md5))
  sprintf("与上一次运行逐文件 md5 比较：%d / %d 相同", same, length(new_md5))
}
if (!is.null(prev_md5)) add_check("复跑一致：R 产物 md5 与上一次运行相同（rng.seed = 1）",
                                  all(unname(prev_md5) == unname(new_md5)), rerun_txt)

## ---------------------------------------------------------------- md
vers <- function(p) as.character(packageVersion(p))
L <- character(0)
A <- function(...) L <<- c(L, ...)
A("# 09_gate1b 任务 2（Q2c）：作者 DAPAR 流程分支", "",
  "由 `02_author_branch.R`（主流程）生成；⑥ 节由 `02_author_compare.py` 填写。只报数。", "")

A("## ① 输入与环境", "",
  sprintf("- 输入：`%s`（%d 行 × %d 列；%d 个 run 全部可由 07 `parse_design` 三个正则解析为 5 tp × 6 FR × 4 rep）。",
          IN_PG_REL, nrow(raw), ncol(raw), length(runs)),
  "- 并列表对照：`analysis/08_gate1/02_events.tsv`（表 E-C，由 `02_author_compare.py` 读取）。",
  sprintf("- 作者脚本：`%s`、`translocation_plots.R`（同目录）。", sub(paste0(REPO, "/"), "", AUTH_D, fixed = TRUE)),
  sprintf("- R：%s；Rscript = `<scratchpad>/envs/rdapar/bin/Rscript`（micromamba 环境，未装新包）。", R.version.string),
  sprintf("- 包版本：limma %s；impute %s；metap %s；Biobase %s；DAPAR %s（未 `library`，见下）；dplyr %s、tidyr %s、stringr %s（仅被 DAPAR `limmaCompleteTest` 交叉核对调用）；MSnbase %s（环境内有，本脚本未用）。",
          vers("limma"), vers("impute"), vers("metap"), vers("Biobase"), vers("DAPAR"),
          vers("dplyr"), vers("tidyr"), vers("stringr"), vers("MSnbase")),
  "- DAPAR 取法：`e <- new.env(); lazyLoad(file.path(.libPaths()[1], \"DAPAR\", \"R\", \"DAPAR\"), envir = e)`，再把 `e` 中闭包的环境设为 `e`（否则找不到 `pkgs.require` 等内部函数）。`library(DAPAR)` 因依赖 DAPARdata（代理 403）不可用。主流程底层直接调 limma / impute / stats；DAPAR 1.38 的 `LOESS`、`findMECBlock`、`reIntroduceMEC`、`getQuantile4Imp`、`limmaCompleteTest` 用于逐 fraction 交叉核对，`make.design.1`、`make.contrast` 用于构造 limma 设计矩阵与对比。",
  "- DAPAR 源码行号引用 `<scratchpad>/dapar_1.38_src.txt`（下称 src:）；作者行号 D: = `DAPAR_script_OSM.R`，T: = `translocation_plots.R`。",
  sprintf("- S0：原始 %d 行 → 去 `Genes` 空 %d 行（D:10）→ 取第一个基因名（D:11-12；含 `;` 的 %d 行）→ 去重复基因名 %d 行（D:13，保留首次出现）→ **%d 行**进入流程。四蛋白对应蛋白组：%s。",
          n_raw, n_drop_empty, n_multi_gene, n_drop_dup, n_s0,
          paste(sprintf("%s %s", GENES4, pg_of4), collapse = "、")),
  sprintf("- 非正值（按缺失）%d 个；非空字符串解析失败 %d 个。", n_nonpos, n_parse_fail), "")

A("## ② 处理顺序（实际使用的函数与参数）", "")
proc <- data.frame(
  step = c("S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7a", "S7b", "S8", "S9",
           "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"),
  author = c("D:9-13", "D:21-28, D:40", "D:41", "D:42", "D:43-58", "D:67", "D:68", "D:69-70", "D:71-73",
             "D:76-82", "D:89-107", "T:16-20", "T:23-25", "T:30-31", "T:33-34", "T:36, T:38-47",
             "T:58-148", "T:150-153", "T:157-171", "T:173", "T:179-184"),
  dapar = c("—", "—", "—", "—（1.38 无 mvFilter）", "—",
            "wrapper.normalizeD src:25-46（LOESS 分支 src:42-43）；LOESS src:4-22（overall src:10-13）",
            "wrapper.impute.KNN src:49-73（逐条件 src:59-67；==0→NA src:68）",
            "findMECBlock src:76-93；reIntroduceMEC src:96-105",
            "getQuantile4Imp src:108-114（1.38 的 wrapper.impute.detQuant src:117-143 为 metacell 版，未用）",
            "limmaCompleteTest src:146-194（OnevsOne src:159-173）；make.design src:197-212；make.design.1 src:215-239；make.contrast src:242-301（contrast=1 src:277-286）；formatLimmaResult src:304-352",
            "—", "—", "—", "—", "—", "—", "—", "—", "—", "—", "—"),
  impl = c(
    "`Genes != \"\"` → `gsub(\"^;\",\"\")` → `gsub(\";.*\",\"\")` → `!duplicated()`；行 ID = 基因名",
    "每 fraction 20 列，按 [2min, 8min, 20min, 90min, CTRL] × Rep1–4 排（作者布局：4 个处理条件块在前，CTRL 在 17:20）",
    "`log2(x)`；≤ 0 按缺失（本数据 0 个）",
    sprintf("旧版 mvFilter 逻辑：保留至少一个条件中非缺失数 ≥ %d 的行（按旧版逻辑实现）", TH),
    "每对比 k：`KEEP` ⇔ CTRL 非缺失 ≥ 3 或 条件 k 非缺失 ≥ 3",
    sprintf("`limma::normalizeCyclicLoess(x, method = \"fast\", span = %s)`（iterations = 3、weights = NULL 为库默认；NA 由 `loessFit` 排除）", SPAN),
    sprintf("逐条件 4 列 `impute::impute.knn(x, k = %d, rowmax = %s, colmax = %s, maxp = %d, rng.seed = %d)`；之后 `x[x == 0] <- NA`（src:68）", KNN_K, ROWMAX, COLMAX, MAXP, RNG_SEED),
    "MEC = 某条件 4 个 rep 在 LOESS 后矩阵中全缺的 (行, 条件) 块；KNN 后把这些块重新置 NA",
    sprintf("每列 `stats::quantile(x, %s, na.rm = TRUE)`（type 7）× %d；旧版 `impute.detQuant(qData, values)`：每列剩余 NA ← 该列值（按旧版逻辑实现）", QVAL, FACTOR),
    "每对比 k：输入 = 条件 k 4 列 + CTRL 4 列，只取 KEEP 行；设计矩阵 = `make.design.1`（条件指示矩阵，非配对）；对比 = `make.contrast(contrast = 1)` → `( ConditionA )/ 1 -( ConditionB )/ 1`（tp − CTRL）；`limma::lmFit → limma::contrasts.fit → limma::eBayes`（默认参数）；`p.adjust(P, \"BH\", n = length(P))` 在该 fraction × 对比的 KEEP 行内",
    "6 个 fraction 的填补后 log2 矩阵按基因合并为 120 列（没过 S3 的块为 NA）；`rowSums(is.na) < 120` 过滤（D:105）",
    "取该 tp 与 CTRL 各 24 列，`2^`（含填补值）",
    "每 fraction 4 个 rep 的算术均值（线性）",
    "份额 = fraction 均值 / 6 个 fraction 均值之和（tp、CTRL 分别）",
    "MS = |份额_tp − 份额_CTRL| 的 6 个值的均值",
    "`match(max(x), x)`；`maxN`（原样）+ `match`",
    "fraction max1 / max2 在该 tp 的 limma 表中的 P_Value（作者第 3 列）；蛋白不在表 → 1",
    "`metap::sumlog(c(p1, p2))$p`",
    "limegreen / dodgerblue / coral / gray，规则原样",
    "`p.adjust(pval_combi, \"BH\")`，在全部可算蛋白上（见 ⑦ (g)）",
    "标注条件 `FDR < 0.05 & MS > 0.1`（严格不等号，原样）"))
colnames(proc) <- c("步", "作者行号", "DAPAR 1.38 源码行号", "本轮实际函数与参数")
A(md_table(proc), "")

A("## ③ 每 fraction 的过滤与填补统计", "",
  sprintf("每 fraction 输入 %d 行 × 20 列。observed = LOESS 后仍为观测值的格；KNN = impute.knn 填补且未被 MEC 重置的格（含 impute.knn 内部退回列均值的格，见 KNN_colmean）；MEC_blocks = 某条件 4 个 rep 全缺的 (行, 条件) 块数；detQuant = 2.5%% 分位数填补的格（= MEC 格 + impute.knn 输出为 0 被 src:68 置 NA 的格）。impute.knn 退回列均值的行（逐条件 5 次调用之和）分两类：rowmax 类 = 缺失比例 > rowmax = 0.99（即 4/4 缺失）的行，impute 发警告并用列均值填，随后 S7a 作为 MEC 重新置 NA、由 detQuant 填；newmiss 类 = KNN 计算后仍为 NA 的行（无警告），用列均值填，这些格留在最终矩阵中、来源记 KNN（KNN_colmean 列为其格数）。", n_s0), "")
t3 <- data.frame(fraction = FRSTAT$fraction, S3_kept = FRSTAT$S3_kept, observed = FRSTAT$observed,
  KNN = FRSTAT$KNN, KNN_colmean = FRSTAT$knn_meanimp_cells_newmiss, MEC_blocks = FRSTAT$MEC_blocks,
  detQuant = FRSTAT$detQuant, colmean_rows_rowmax = FRSTAT$knn_meanimp_rows_rowmax,
  colmean_rows_newmiss = FRSTAT$knn_meanimp_rows_newmiss,
  check = sprintf("%d+%d+%d = %d = %d×20 %s", FRSTAT$observed, FRSTAT$KNN, FRSTAT$detQuant,
                  FRSTAT$observed + FRSTAT$KNN + FRSTAT$detQuant, FRSTAT$S3_kept, ifelse(FRSTAT$sum_ok, "✓", "✗")))
tot <- data.frame(fraction = "合计", S3_kept = sum(FRSTAT$S3_kept), observed = sum(FRSTAT$observed),
  KNN = sum(FRSTAT$KNN), KNN_colmean = sum(FRSTAT$knn_meanimp_cells_newmiss), MEC_blocks = sum(FRSTAT$MEC_blocks),
  detQuant = sum(FRSTAT$detQuant), colmean_rows_rowmax = sum(FRSTAT$knn_meanimp_rows_rowmax),
  colmean_rows_newmiss = sum(FRSTAT$knn_meanimp_rows_newmiss), check = "")
A(md_table(rbind(t3, tot)), "")
A(sprintf("- rowmax 类退回列均值的行数 = MEC 块数（每 fraction 相等：%s），且等于 impute.knn 警告 \"N rows with more than 99 %% entries missing; mean imputation used for these rows\" 中 N 之和。",
          if (all(FRSTAT$knn_meanimp_rows_rowmax == FRSTAT$MEC_blocks)) "是" else "否"),
  sprintf("- detQuant 格 = 4 × MEC 块 + impute.knn 输出为 0 的格：%s。impute.knn 输出为 0 的格每 fraction 为 %s（只可能来自 knnimp.split 切出的 ≤ k 行小簇：库代码对其调 `meanimp(x[index, ])` 时缺失格已置 0、识别不到缺失；复刻计数小簇内缺失格 %s），故 src:68 的 `== 0 → NA` 在本数据上无实际作用。",
          paste(sprintf("%s %d = 4×%d + %d", FRSTAT$fraction, FRSTAT$detQuant, FRSTAT$MEC_blocks, FRSTAT$detQuant_zero), collapse = "；"),
          paste(FRSTAT$knn_zero_cells, collapse = " / "), paste(FRSTAT$knn_smallcluster_cells, collapse = " / ")),
  sprintf("- S8 limma 行数（KEEP 行）：%s。", paste(sprintf("%s 2min %d / 8min %d / 20min %d / 90min %d", FRSTAT$fraction,
          FRSTAT$KEEP_2min, FRSTAT$KEEP_8min, FRSTAT$KEEP_20min, FRSTAT$KEEP_90min), collapse = "；")),
  sprintf("- detQuant 值（每列 2.5%% 分位数 × 1，log2）范围：%s。120 个值见 `02_author_detquant_values.tsv`。",
          paste(vapply(FRS, function(fr) { v <- DETQ$shiftedImpVal_log2[DETQ$fraction == fr]
            sprintf("%s %.3f–%.3f", fr, min(v), max(v)) }, ""), collapse = "；")),
  sprintf("- S9：合并表（D:105 过滤后）%d 行 × 120 列（%d 行在 6 个 fraction 都没过 S3 而被去掉）；其中 120 列全非 NA（6 个 fraction 都过 S3）的 %d 行进入 T1–T10。",
          n_merged, n_s0 - n_merged, n_complete),
  sprintf("- S3 保留次数分布（9607 行中，在几个 fraction 过 S3）：%s。",
          paste(sprintf("%d 个 %d 行", as.integer(names(n_pass_tab)), as.integer(n_pass_tab)), collapse = "；")), "")
A("逐条件明细：", "")
cs <- COND_STAT[, c("fraction", "condition", "missing_cells", "MEC_blocks", "knn_meanimp_rows_rowmax",
                    "knn_warning_rows", "knn_meanimp_rows_newmiss", "knn_meanimp_cells_newmiss",
                    "knn_smallcluster_cells", "knn_cluster_lines")]
A(md_table(cs), "")

A("## ④ 四蛋白每 fraction 的状态", "",
  "来源代码：O = observed，K = KNN（M = KNN 步骤中 impute.knn 退回列均值的格，tsv 中来源仍记 KNN），D = detQuant，- = filtered_out（该 fraction 没过 S3）；每格 4 个字符依次为 Rep1–Rep4。`n_obs` 为该 fraction 各条件原始非缺失数（CTRL/2/8/20/90）；过 S3 需至少一个条件 ≥ 3。KEEP 列为该蛋白是否进入该 fraction 对应对比（2/8/20/90min vs CTRL）的 limma 表。", "")
code_of <- c(observed = "O", KNN = "K", detQuant = "D", filtered_out = "-")
code_cell <- function(g, r, fr) {
  s0c <- Sall[g, r]
  if (s0c == "KNN" && isTRUE(COLMEAN[[fr]][g, r])) "M" else code_of[[s0c]]
}
for (g in GENES4) {
  n_cm <- sum(vapply(FRS, function(fr) if (g %in% rownames(COLMEAN[[fr]])) sum(COLMEAN[[fr]][g, ]) else 0L, 0))
  fail_fr <- FRS[!vapply(FRS, function(fr) g %in% rownames(FINAL[[fr]]), TRUE)]
  A(sprintf("**%s**（%s）：原始非空格 %d；observed %d、KNN %d（其中列均值退回 %d）、detQuant %d、filtered_out %d（其中原始非空 %d）。没过 S3 的 fraction：%s。",
            g, pg_of4[[g]], raw_nonempty4[[g]], sum(cells$gene == g & cells$source == "observed"),
            sum(cells$gene == g & cells$source == "KNN"), n_cm, sum(cells$gene == g & cells$source == "detQuant"),
            sum(cells$gene == g & cells$source == "filtered_out"), g_fo_raw[[g]],
            if (length(fail_fr)) paste(fail_fr, collapse = "、") else "无"), "")
  rows <- data.frame()
  for (fr in FRS) {
    passed <- g %in% rownames(FINAL[[fr]])
    nob <- vapply(c("CTRL", TPS), function(cn) sum(!is.na(L2[g, vapply(REPS, function(rp) run_of(cn, fr, rp), "")])), 0)
    codes <- vapply(c("CTRL", TPS), function(cn) paste(vapply(REPS, function(rp) code_cell(g, run_of(cn, fr, rp), fr), ""), collapse = ""), "")
    keepk <- if (passed) paste(ifelse(MVK[[fr]]$KEEP[g, ], "Y", "N"), collapse = "/") else "—"
    rows <- rbind(rows, data.frame(fraction = fr, S3 = if (passed) "过" else "没过",
      n_obs = paste(nob, collapse = "/"), CTRL = codes[["CTRL"]], `2min` = codes[["2min"]],
      `8min` = codes[["8min"]], `20min` = codes[["20min"]], `90min` = codes[["90min"]],
      `KEEP 2/8/20/90` = keepk, check.names = FALSE))
  }
  A(md_table(rows), "")
}

A("## ⑤ 结果表（`02_author_branch.tsv`）", "",
  sprintf("T1–T10 的蛋白集合 = 合并表中 6 个 fraction 都过 S3 的 %d 行（BH 的 n）。没过 S3 的蛋白（某 fraction 为 NA）在作者代码中无法计算（见 ⑦ (g)），本表留空。MS、份额 4 位小数；p 3 位有效数字；全精度值见 tsv。",
          n_complete),
  sprintf("四蛋白是否在这 %d 行内：%s。", n_complete, paste(vapply(GENES4, function(g) {
    ff <- FRS[!vapply(FRS, function(fr) g %in% rownames(FINAL[[fr]]), TRUE)]
    if (length(ff)) sprintf("%s 否（%s 没过 S3）", g, paste(ff, collapse = "、")) else sprintf("%s 是", g) }, ""), collapse = "；")), "")
b5 <- data.frame(gene = br$gene, timepoint = br$timepoint,
  MS = f4(suppressWarnings(as.numeric(br$MS_author_pipeline))), max1 = br$MS_max1, max2 = br$MS_max2,
  share_max1_tp = f4(suppressWarnings(as.numeric(br$share_max1_tp))),
  share_max1_ctrl = f4(suppressWarnings(as.numeric(br$share_max1_ctrl))),
  dir = br$direction_max1, p_max1 = g3(suppressWarnings(as.numeric(br$p_max1))),
  p_max2 = g3(suppressWarnings(as.numeric(br$p_max2))), pval_combi = g3(suppressWarnings(as.numeric(br$pval_combi))),
  FDR = g3(suppressWarnings(as.numeric(br$pval_combi_FDR))), color = br$color_category,
  `MS>0.1` = br$MS_gt_0.1, `FDR<0.05` = br$FDR_lt_0.05, n_FDR = br$n_proteins_in_FDR, check.names = FALSE)
A(md_table(b5), "")
A("全蛋白计数（每个 tp）：", "")
A(md_table(PROT_SUM), "",
  "（n_label_T181 = 满足作者 T:181 标注条件 `FDR < 0.05 & MS > 0.1` 的蛋白数；n_p1_eq_1 / n_p2_eq_1 = max1 / max2 fraction 的 limma 表中没有该蛋白而取 p = 1 的蛋白数。）", "")

A("## ⑥ 并列表（`02_author_compare.tsv`）", "",
  "<!-- SECTION6_BEGIN -->", "（由 `02_author_compare.py` 填写；尚未运行。）", "<!-- SECTION6_END -->", "")

A("## ⑦ 与作者实现的差别清单", "",
  "- (a) 输入：本轮是 120 个单独 DIA-NN 搜库的 pg_matrix 合并表（`proteome_pg_matrix_120.tsv`，线性 PG 定量），不是作者一次 Spectronaut 合并搜库的蛋白报告（D:8-9）。S0 的去空、取第一个基因名、去重按 D:10-13 照做（本数据无以 `;` 开头的基因名）。",
  "- (b) 条件：EGF 五个时间点（CTRL, 2min, 8min, 20min, 90min）；对比 4 个（2/8/20/90min vs CTRL）。作者 OSM 脚本为 5 个渗透压条件；translocation 作者只做 2min vs CTRL，本轮对四个 EGF 时间点各做一遍（把 T:16 的 `\"2min\"` 换成各 tp，其余代码不变）。",
  sprintf("- (c) KNN 随机种子：作者 `rng.seed = sample(seq_len(1000), 1)`（src:65，随机，未设种子）；本轮每次调用固定 `rng.seed = %d`。每个条件的行数 > maxp = 1500，impute 走二分聚类（用随机起点），故种子影响 KNN 结果；作者的结果本身不可复现。", RNG_SEED),
  "- (d) `mvFilter`、`impute.detQuant`：DAPAR 1.38 已无（改为 metacell 机制），按旧版逻辑实现：mvFilter(atLeastOneCond, th = 3) = 至少一个条件非缺失数 ≥ 3 保留；impute.detQuant(qData, values) = 每列剩余 NA 用该列的 shiftedImpVal 替换。",
  sprintf("- (e) LOESS：以 DAPAR 1.38 源码为准（src:10-13：`limma::normalizeCyclicLoess(x = qData, method = \"fast\", span = span)`，span 默认 0.7），配 limma %s。作者 2020 年所用 DAPAR / limma 版本的 LOESS 实现是否完全相同无法核对。", vers("limma")),
  "- (f) T1 的 `2^` 输入是填补后的 log2 值（含 KNN 与 detQuant 填补值），与作者相同（作者 T:7 读的 merged_table 即 D:93-107 的填补后矩阵）。",
  sprintf("- (g) 【方案未定义处】没过 S3 的 (蛋白, fraction) 块在合并表中为 NA（D:93 `merge(..., all = T)`）。作者 T:24-34 对含 NA 的行得到全 NA 的份额与 MS，T:47 的 `maxN` 随即报错（`sort.int: index 5 outside bounds`，本环境实测），即作者代码在含此类行的表上无法运行完。本轮把 T1–T10 限于 120 列全非 NA 的 %d 行（BH 的 n = %d），其余蛋白的 MS / p / FDR 留空。四蛋白中受影响的见 ④、⑤。", n_complete, n_complete),
  "- (h) T:56 作者用 `full_table$Gene.names` 作行名；OSM 版 D:107 写出的 merged_table 不含 `Gene.names` 列（基因名在行名中），作者 EGF 版 DAPAR 脚本不在仓库。本轮按方案 T6 用基因名（合并表行 ID）匹配 limma 表的行名（作者 csv 的 `X` 列）。",
  "- (i) DAPAR 1.38 的 `wrapper.impute.KNN` 另有按 metacell 把 \"Missing MEC\" 置 NA 并更新 metacell（src:69-71）；本轮不建 metacell，用作者显式的 findMECBlock / reIntroduceMEC（D:69-70）实现同一效果。src:68 的 `== 0 → NA` 照 1.38 保留（本数据上 0 格，无实际作用，见 ③）；作者所用旧版 wrapper.impute.KNN 是否含此行无法离线核对。",
  "- (j) limma：主流程直接调 `limma::lmFit → contrasts.fit → eBayes`，设计矩阵与对比字符串用 DAPAR 1.38 的 `make.design.1` / `make.contrast` 生成；与 DAPAR 1.38 `limmaCompleteTest` 整体调用的结果逐 fraction × 对比核对一致（见自检）。BH 照 D:81 在该 fraction × 对比的 KEEP 行内做。",
  "- (k) 输出形式：limma 表按方案写成 tsv（列 gene, logFC, P_Value, BH），对应作者 D:82 的 csv（行名 = 基因名、logFC、P_Value、BH_pvalue）；T6 取的 \"第 3 列\" = P_Value（未校正）。", "")

A("## ⑧ 跳过项", "",
  "- D:83-85 的 .rnk 文件、D:89-91 / D:97-99 的每 fraction 条件均值 / 2^ / SD 表、D:101-102 的 `cbind(mv_keep, prot_imp_2)` 表：不影响 T1–T10，未写出。",
  "- T:179-184 的 ggplot 火山图：未画；标注条件 `FDR < 0.05 & MS > 0.1` 以计数与 ⑤ 的 `MS_gt_0.1`、`FDR_lt_0.05` 列给出。",
  "- `设计v1_final_2026-10-06.md`：缺，按 00_design 任务书执行。",
  "- 无其它跳过；无新装 R 包。", "")

A("## 自检（R 部分）", "")
A(md_table(SELF), "", sprintf("- 复跑：%s。", rerun_txt), "")
A("交叉核对明细（每 fraction）：", "")
xk <- XCHK; xk$limma_vs_limmaCompleteTest_maxdiff <- sprintf("%.1e", xk$limma_vs_limmaCompleteTest_maxdiff)
xk$observed_vs_loess_maxdiff <- sprintf("%.1e", xk$observed_vs_loess_maxdiff)
A(md_table(xk), "")

## 若 md 已存在且 ⑥ 已由 Python 填写，保留 ⑥ 的内容
if (file.exists(F_MD)) {
  old <- readLines(F_MD, warn = FALSE, encoding = "UTF-8")
  b <- which(old == "<!-- SECTION6_BEGIN -->"); e2 <- which(old == "<!-- SECTION6_END -->")
  nb <- which(L == "<!-- SECTION6_BEGIN -->"); ne <- which(L == "<!-- SECTION6_END -->")
  if (length(b) == 1 && length(e2) == 1 && e2 > b) L <- c(L[seq_len(nb)], old[(b + 1):(e2 - 1)], L[ne:length(L)])
}
writeLines(L, F_MD, useBytes = TRUE)

cat("\n=== SELF-CHECK (R) ===\n"); print(SELF, right = FALSE)
cat("\n=== FRSTAT ===\n"); print(FRSTAT)
cat("\n=== branch ===\n"); print(br)
cat("\n", rerun_txt, "\n")
if (any(SELF$result != "PASS")) quit(status = 1)
