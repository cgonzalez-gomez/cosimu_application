# devtools::install_git("git@gitlab.signiatx.eu:signia/bioinfo/packages/coinf.git",
#                       git = "external",upgrade = F,force = T,dependencies = T)
# rstudioapi::restartSession()
library(coinf)
library(tictoc)
get_list_deg <- function(fc_vect, abs_thresh = 1.5, 
                         pval = NA,use.pval = FALSE, thresh_pval = 0.05, 
                         use.top=FALSE,topN = NA){
  idx <- names(fc_vect)
  if(is.null(idx)){
    stop("fc_vect needs to be a named vector")
  }
  if(use.top){
    if(topN < length(fc_vect)/2){
      up <- idx[order(fc_vect,na.last = T,decreasing = T)[1:topN]]
      down <- idx[order(fc_vect,na.last = T,decreasing = F)[1:topN]]
    }else{
      stop("All the genes are selected (top N is too big), use a lfc threshold of 0
           in this case")
    }
  }else if(use.pval){
    up <- idx[fc_vect >= abs_thresh & !is.na(fc_vect) & pval < thresh_pval & !is.na(pval)]
    down <- idx[fc_vect <= -abs_thresh & !is.na(fc_vect) & pval < thresh_pval & !is.na(pval)]
  }else{
    up <- idx[fc_vect >= abs_thresh & !is.na(fc_vect)]
    down <- idx[fc_vect <= -abs_thresh & !is.na(fc_vect)]
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
if(!is.na(args[5])){
  mat_pval <- readRDS(args[5])
}

paired_est <- as.logical(args[6])
if(paired_est){
  params <- readRDS(args[7])
}
permuted_pval = as.logical(args[8])
topN <- NA
if(!is.na(args[9])){
  topN <- as.numeric(args[9])
}
gene_names <- rownames(mat_profiles)
pert_names <- colnames(mat_profiles)
# Filter the query with the topN not the DB
# L3 <- coinf::matrixToRankedList(mat_profiles, topN = topN,ncpus = ncpus) 
L3 <- coinf::matrixToRankedList(mat_profiles, expression_val = T,ncpus = ncpus)

if(!paired_est){
  if(is.na(topN)){
    if(is.null(mat_pval)){
      set.seed(100)
      tic("XSum")
      res_xsum<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        names(fc) <- gene_names
        S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = NA, use.pval = FALSE,
                               thresh_pval = 0.05)
        xsum <- XSum_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                           pert_names = pert_names,list_ranked_R = L3,
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        return(xsum)
      })
      toc()
    }else{
      set.seed(100)
      tic("XSum")
      res_xsum<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        pval <- mat_pval[,i]
        names(fc) <- gene_names
        S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = pval, use.pval = TRUE, thresh_pval = 0.05)
        xsum <- XSum_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                           pert_names = pert_names,list_ranked_R = L3,
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        return(xsum)
      })
      toc()
    }
  }else{
    set.seed(100)
    tic("XSum")
    res_xsum<- purrr::map(seq(1,ncol(mat_queries)),function(i){
      fc <- mat_queries[,i]
      names(fc) <- gene_names
      S_list <- get_list_deg(fc, use.top = TRUE,topN = topN)
      xsum <- XSum_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                         pert_names = pert_names,list_ranked_R = L3,
                         ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
      return(xsum)
    })
    toc()
  }
}else{
  if(is.na(topN)){
    if(is.null(mat_pval)){
      set.seed(100)
      tic("XSum")
      res_xsum<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        names(fc) <- gene_names
        sel <- which(params$ind_base==i)
        S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = NA, use.pval = FALSE, thresh_pval = 0.05)
        xsum <- XSum_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                           pert_names = pert_names[sel],list_ranked_R = L3[sel],
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        return(xsum)
      })
      toc()
    }else{
      set.seed(100)
      tic("XSum")
      res_xsum<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        pval <- mat_pval[,i]
        names(fc) <- gene_names
        sel <- which(params$ind_base==i)
        S_list <- get_list_deg(fc, abs_thresh = 1.5, pval = pval, use.pval = TRUE, thresh_pval = 0.05)
        xsum <- XSum_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                           pert_names = pert_names[sel],list_ranked_R = L3[sel],
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        return(xsum)
      })
      toc()
    }
  }else{
    set.seed(100)
    tic("XSum")
    res_xsum<- purrr::map(seq(1,ncol(mat_queries)),function(i){
      fc <- mat_queries[,i]
      names(fc) <- gene_names
      sel <- which(params$ind_base==i)
      S_list <- get_list_deg(fc, use.top = TRUE,topN = topN)
      xsum <- XSum_score(S_up = S_list$UP, S_down = S_list$DOWN, gene_names = gene_names,
                         pert_names = pert_names[sel],list_ranked_R = L3[sel],
                         ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
      return(xsum)
    })
    toc()
  }
}
saveRDS(res_xsum,paste0(save_path,"_",topN,".rds"))
