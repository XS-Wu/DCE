###########################################################
# Copyright © 2025 Xinsheng Wu, Fudan University
###########################################################
library(idefix)

options(stringsAsFactors = FALSE)
set.seed(20250918)

lvl_names <- list(
  `PrEP Modalities`      = c("每天吃一片药","2-1-1按需吃药","2月注射一次","6月注射一次","皮下植入型药物（持续约4月）"),
  `Cost (per month)`     = c("完全免费，政府全额补贴","政府补贴65%","政府补贴30%","全额自费"),
  `Income (per month)`   = c("远低于平均水平","低于平均水平","大致与平均水平相当","高于平均水平","远高于平均水平"),
  `Adherence Management` = c("自行管理（没有额外提醒）","手机应用程序提醒","短信提醒","匿名互助群"),
  `Accessibility`        = c("医院或诊所","社区卫生服务中心","药房","专为MSM服务的社区组织","快递送药"),
  `Prescription`         = c("不需要等化验结果，当天就能拿到处方","等化验结果确认符合用药条件后再拿处方",
                             "需要做相同的化验检查并上传结果，通过线上平台获取处方"),
  `Side effects`         = c("轻微（如恶心、食欲不振等胃肠道症状）","中等（如头痛、腹泻、无力）","严重（如骨密度降低、肾功能受损）")
)

atts   <- names(lvl_names)
lvls   <- vapply(lvl_names, length, integer(1))      
coding <- rep("D", length(lvls))                     
P      <- sum(lvls - 1L)                            
NSETS  <- 60
NALTS  <- 2


prior_list <- list(
  `PrEP Modalities`      = c(+0.10, +0.20, +0.40, +0.25),    
  `Cost (per month)`     = c(-1.935, -1.935, -4.683),        
  `Income (per month)`   = c(+0.80, +1.60, +2.50, +4.00),    
  `Adherence Management` = c(-0.149, -0.253, -0.014),        
  `Accessibility`        = c(-0.144, -0.031, -0.170, +0.120),
  `Prescription`         = c(-0.25, -0.50),                  
  `Side effects`         = c(-1.00, -5.00)                   
)

beta0 <- unlist(lapply(atts, function(a) {
  v <- prior_list[[a]]; stopifnot(length(v) == lvls[[a]] - 1L); v
}))
stopifnot(length(beta0) == P)

`%||%` <- function(x, y) if (is.null(x)) y else x

make_blocks_base <- function(ques_df, n_blocks = 3, seed = 2025) {
  stopifnot(n_blocks >= 1)
  ids <- sort(unique(ques_df$obsID))
  stopifnot(length(ids) %% n_blocks == 0)
  set.seed(seed)
  ids <- sample(ids)                              
  split(ids, rep(1:n_blocks, each = length(ids)/n_blocks))
}

run_one_design <- function(seed, label,
                           n_sets = 60, n_alts = 2,
                           n_start = 8,  max_iter = 120,
                           do_blocks = TRUE, n_blocks = 3) {
  set.seed(seed); flush.console()
  flush.console()
  
  atts_raw  <- names(lvl_names)
  lvls_raw  <- vapply(lvl_names, length, integer(1))
  coding    <- rep("D", length(lvls_raw))
  P         <- sum(lvls_raw - 1L)
  att_safe  <- sprintf("X%02d", seq_along(lvls_raw))
  lvls_safe <- lvls_raw; names(lvls_safe) <- att_safe
  
  beta0 <- unlist(lapply(atts_raw, function(a) {
    v <- switch(a,
                "PrEP Modalities"       = c(+0.10, +0.20, +0.30, +0.25),
                "Cost (per month)"      = c(-0.40, -0.80, -1.20),
                "Income (per month)"    = c(+0.05, +0.10, +0.15, +0.20),
                "Adherence Management"  = c(+0.15, +0.12, +0.08),
                "Accessibility"         = c(+0.10, +0.12, +0.20, +0.25),
                "Prescription"          = c(-0.15, -0.10),
                "Side effects"          = c(-0.50, -1.10)
    ); v
  }))
  stopifnot(length(beta0) == P)
  
  flush.console()
  des_obj <- try(idefix::CEA(
    lvls      = unname(lvls_safe),
    coding    = coding,
    n.sets    = n_sets,
    n.alts    = n_alts,
    par.draws = beta0,
    optim     = "D",
    parallel  = FALSE,
    n.start   = n_start,
    max.iter  = max_iter
  ), silent = TRUE)
  if (inherits(des_obj, "try-error")) stop("CEA 失败（安全模式下不回退 Modfed）。")
  flush.console()
  
  get_design_mat <- function(x) {
    if (is.list(x) && "BestDesign" %in% names(x)) return(x$BestDesign$design)
    if (is.list(x) && "design"     %in% names(x)) return(x$design)
    if (is.matrix(x)) return(x)
    stop("无法从对象中提取设计矩阵。")
  }
  des_mat <- get_design_mat(des_obj)
  
  flush.console()
  eval <- idefix::EvaluateDesign(des = des_mat, par.draws = matrix(beta0,1), n.alts = n_alts)
  flush.console()
  
  flush.console()
  dec <- idefix::Decode(des = des_mat, n.alts = n_alts, lvl.names = lvl_names, coding = coding)
  char_tab <- if (is.list(dec)) (dec$des %||% dec$design %||% dec$char.des %||% dec[[1]]) else dec
  n_sets_eff <- nrow(char_tab) / n_alts
  ques <- data.frame(
    obsID = rep(seq_len(n_sets_eff), each = n_alts),
    alt   = rep(1:n_alts, times = n_sets_eff),
    char_tab, check.names = FALSE
  )
  flush.console()
  
  if (isTRUE(do_blocks)) {
    flush.console()
    blks <- make_blocks_base(ques, n_blocks = n_blocks, seed = seed + 1L)
    for (b in seq_along(blks)) {
      ids <- blks[[b]]
      sub <- ques[ques$obsID %in% ids, ]
      ord <- sample(ids)
      sub <- sub[order(match(sub$obsID, ord), sub$alt), ]
      fn  <- sprintf("Blocks_%s_Block%02d.csv", label, b)
      write.csv(sub, fn, row.names = FALSE, fileEncoding = "UTF-8")
    }
    flush.console()
  }
  
  flush.console()
  fn_q <- sprintf("DCE_60x2_idefix_questionnaire_%s.csv", label)
  fn_x <- sprintf("DCE_60x2_idefix_designMatrix_%s.csv",   label)
  write.csv(ques,   fn_q, row.names = FALSE, fileEncoding = "GBK")
  write.csv(des_mat, fn_x, row.names = FALSE)   
  flush.console()
  
  invisible(list(des_mat = des_mat, ques = ques, eval = eval,
                 files = list(questionnaire = fn_q, design = fn_x)))
}


A <- run_one_design(seed = 20250918, label = "A",
                    n_start = 12, max_iter = 200,
                    do_blocks = TRUE, n_blocks = 3)
B <- run_one_design(seed = 20251999, label = "B",
                    n_start = 12, max_iter = 200,
                    do_blocks = TRUE, n_blocks = 3)

attr_tbl <- do.call(rbind, lapply(names(lvl_names), function(a){
  data.frame(Attribute = a, Level = lvl_names[[a]], check.names = FALSE)
}))

as_overlap_hist <- function(ov) {
  vals <- tryCatch({
    if (is.list(ov)) {
      unlist(ov, recursive = TRUE, use.names = FALSE)
    } else if (is.matrix(ov) || is.data.frame(ov)) {
      as.vector(as.matrix(ov))
    } else {
      as.vector(ov)
    }
  }, error = function(e) numeric(0))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return("<NA>")
  vals <- as.integer(round(vals))
  tb <- table(vals)
  paste(paste0(names(tb), ":", as.integer(tb)), collapse = " | ")
}

mk_diag_row <- function(x, lab) {
  ev <- x$eval
  dberr <- as.numeric(ev$DB.error)[1]       
  ortho <- as.numeric(ev$Orthogonality)[1]
  ovh   <- as_overlap_hist(ev$level.overlap)
  n_tasks <- nrow(x$ques) / length(unique(x$ques$alt))
  n_alts  <- length(unique(x$ques$alt))
  P <- sum(vapply(lvl_names, function(v) length(v) - 1L, integer(1)))
  data.frame(
    Label           = lab,
    N_tasks         = n_tasks,
    N_alts_per_task = n_alts,
    N_params        = P,
    DB_error_local  = round(dberr, 6),
    Orthogonality   = round(ortho, 4),
    Overlap_hist    = ovh,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}



write_attr_levels <- function(lvl_names, fn="Table_Attributes_and_Levels.csv"){
  rows <- do.call(rbind, lapply(names(lvl_names), function(a){
    data.frame(Attribute=a, Level=lvl_names[[a]], check.names=FALSE)
  }))
  write.csv(rows, fn, row.names=FALSE, fileEncoding="UTF-8"); rows
}

write_priors <- function(lvl_names, beta0, fn="Table_Priors.csv"){
  params <- unlist(lapply(names(lvl_names), function(a){
    L <- length(lvl_names[[a]]); if (L<=1) return(character(0))
    paste0(a, "_", 2:L)
  }))
  stopifnot(length(params)==length(beta0))
  tb <- data.frame(Param=params, Beta0=as.numeric(beta0), check.names=FALSE)
  write.csv(tb, fn, row.names=FALSE, fileEncoding="UTF-8"); tb
}

summarize_design <- function(des_csv, n_alts=2, beta=beta0){
  M <- as.matrix(read.csv(des_csv, check.names=FALSE))
  stopifnot(nrow(M) %% n_alts == 0)
  T <- nrow(M)/n_alts
  DX <- M[seq(2, by=2, length.out=T), ] - M[seq(1, by=2, length.out=T), ]
  eta <- as.numeric(DX %*% beta)
  p <- 1/(1+exp(-eta)); w <- p*(1-p); w[w<1e-10] <- 1e-10
  Z <- DX * sqrt(w); I <- crossprod(Z) + diag(1e-10, ncol(Z))
  d_error <- -as.numeric(determinant(I, TRUE)$modulus)/ncol(I)
  DB <- exp(d_error)
  z <- sweep(DX, 2, beta, `*`)
  dom_BoverA <- sum(apply(z >= -1e-12, 1, all) & apply(z > 1e-12, 1, any))
  dom_AoverB <- sum(apply(z <=  1e-12, 1, all) & apply(z < -1e-12, 1, any))
  dup_pairs  <- sum(apply(DX==0, 1, all))
  list(DB.error=DB, d.error=d_error,
       Dominance_B_over_A=dom_BoverA, Dominance_A_over_B=dom_AoverB,
       Duplicate_pairs=dup_pairs)
}

prior_sensitivity <- function(des_csv, scales=c(0.5,1,1.5)){
  M <- as.matrix(read.csv(des_csv, check.names=FALSE))
  T <- nrow(M)/2; DX <- M[seq(2,by=2,length.out=T),] - M[seq(1,by=2,length.out=T),]
  f <- function(b){
    eta <- as.numeric(DX %*% b); p <- 1/(1+exp(-eta)); w <- p*(1-p); w[w<1e-10] <- 1e-10
    Z <- DX * sqrt(w); I <- crossprod(Z) + diag(1e-10, ncol(Z))
    d <- -as.numeric(determinant(I, TRUE)$modulus)/ncol(I); exp(d)
  }
  data.frame(Scale=scales, DB.error=sapply(scales, function(s) f(s*beta0)))
}





