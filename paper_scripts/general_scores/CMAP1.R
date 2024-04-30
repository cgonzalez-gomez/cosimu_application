# devtools::install_git("git@gitlab.signiatx.eu:signia/bioinfo/packages/coinf.git",
#                       git = "external",upgrade = F,force = T,dependencies = T)
# rstudioapi::restartSession()
library(coinf)
library(tictoc)
get_list_deg <- function(fc_vect, abs_thresh = 1.5, pval = NA,use.pval = FALSE, thresh_pval = 0.05){
  idx <- names(fc_vect)
  if(is.null(idx)){
    stop("fc_vect needs to be a named vector")
  }
  if(!use.pval){
    up <- idx[fc_vect >= abs_thresh & !is.na(fc_vect)]
    down <- idx[fc_vect <= -abs_thresh & !is.na(fc_vect)]
  }else{
    up <- idx[fc_vect >= abs_thresh & pval < thresh_pval & !is.na(fc_vect) & !is.na(pval)]
    down <- idx[fc_vect <= -abs_thresh & pval < thresh_pval & !is.na(fc_vect) & !is.na(pval)]
  }
  return(list("UP"=up[!is.na(up)], "DOWN"=down[!is.na(down)]))
}
args = commandArgs(trailingOnly=TRUE)
args[args==""] <- NA
mat_queries <- readRDS(args[1]) # independent
mat_profiles <- readRDS(args[2]) # dependent
ncpus <- as.numeric(args[3]) #ncpus
save_path <- args[4]
if(is.na(ncpus)) ncpus <- 1

mat_pval <- NULL
if(!is.null(args[5])){
  mat_pval <- readRDS(args[5])
}

paired_est <- as.logical(args[6])
if(paired_est){
  params <- readRDS(args[7])
}
permuted_pval = as.logical(args[8])
# topN <- as.numeric(args[9])
# if(is.na(topN)){
#   topN <- nrow(mat_profiles)
# }
gene_names <- rownames(mat_profiles)
pert_names <- colnames(mat_profiles)
L1 <- coinf::matrixToRankedList(mat_profiles,ncpus = ncpus)
if(!paired_est){
if(is.null(mat_pval)){
  set.seed(100)
  tic("CMAP1")
  res_cmap1 <- purrr::map(seq(1,ncol(mat_queries)),function(i){
    fc <- mat_queries[,i]
    names(fc) <- gene_names
    S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = NA, use.pval = FALSE, thresh_pval = 0.05)
    c1 <- cmap1_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                      pert_names = pert_names,list_ranked_R = L1,ncpus = ncpus,
                      return_raw = F, normalize = T,permuted_pval = permuted_pval,permute_nb = 100)
    return(c1)
  })
}else{
  set.seed(100)
  tic("CMAP1")
  res_cmap1 <- purrr::map(seq(1,ncol(mat_queries)),function(i){
    fc <- mat_queries[,i]
    names(fc) <- gene_names
    pval <- mat_pval[,i]
    names(fc) <- gene_names
    S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = pval, use.pval = TRUE, thresh_pval = 0.05)
    c1 <- cmap1_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                      pert_names = pert_names,list_ranked_R = L1,ncpus = ncpus,
                      return_raw = F, normalize = T,permuted_pval = permuted_pval,permute_nb = 100) 
    return(c1)
  })
  toc()
}
}else{
if(is.null(mat_pval)){
  set.seed(100)
  tic("CMAP1")
  res_cmap1 <- purrr::map(seq(1,ncol(mat_queries)),function(i){
    fc <- mat_queries[,i]
    names(fc) <- gene_names
    sel <- which(params$ind_base==i)
    S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = NA, use.pval = FALSE, thresh_pval = 0.05)
    c1 <- cmap1_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                      pert_names = pert_names[sel],list_ranked_R = L1[sel],ncpus = ncpus,
                      return_raw = F, normalize = T,permuted_pval = permuted_pval,permute_nb = 100)
    return(c1)
  })
}else{
  set.seed(100)
  tic("CMAP1")
  res_cmap1 <- purrr::map(seq(1,ncol(mat_queries)),function(i){
    fc <- mat_queries[,i]
    names(fc) <- gene_names
    pval <- mat_pval[,i]
    sel <- which(params$ind_base==i)
    names(fc) <- gene_names
    S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = pval, use.pval = TRUE, thresh_pval = 0.05)
    c1 <- cmap1_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                      pert_names = pert_names[sel],list_ranked_R = L1[sel],ncpus = ncpus,
                      return_raw = F, normalize = T,permuted_pval = permuted_pval,permute_nb = 100)
    return(c1)
  })
  toc()
}
}
saveRDS(res_cmap1,paste0(save_path,".rds"))
