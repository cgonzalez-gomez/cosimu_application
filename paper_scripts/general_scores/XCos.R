# devtools::install_git("git@gitlab.signiatx.eu:signia/bioinfo/packages/coinf.git",
#                       git = "external",upgrade = F,force = T,dependencies = T)
# rstudioapi::restartSession()
library(coinf)
library(tictoc)

get_deg <- function(fc_vect, abs_thresh = 1.5, pval = NA,use.pval = FALSE, thresh_pval = 0.05, use.top=FALSE,
                    topN=NA){
  if(is.null(names(fc_vect))) stop("fc_vect needs to be a named vector")
  if(use.top){
    if(topN < length(fc_vect)/2){
      fc_vect <- fc_vect[!is.na(fc_vect)]
      fc_vect <- fc_vect[order(fc_vect, decreasing = T)]
      S <- c(utils::head(fc_vect, n = topN), utils::tail(fc_vect, n = topN))
      S <- S[!is.na(S)]
    }else{
      S <- S[!is.na(fc_vect)]
    }
  }else if (use.pval){
    S <- fc_vect[((fc_vect >= abs_thresh & pval < thresh_pval) | (fc_vect <= -abs_thresh & pval < thresh_pval)) & !is.na(fc_vect) & !is.na(pval)]
  }else{
    S <- fc_vect[(fc_vect >= abs_thresh | fc_vect <= -abs_thresh) & !is.na(fc_vect)]
  }
  return(S)
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

paired_est <- as.logical(args[6]) # boolean do we estimate the connectivi VS all the profile or only paired
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
      tic("Xcos")
      res_xcos<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        names(fc) <- gene_names
        S <-get_deg(fc, abs_thresh = 1.5, pval = NA, use.pval = FALSE, thresh_pval = 0.05)
        xcos <- XCos_score(S = S, gene_names = gene_names,
                           pert_names = pert_names,list_ranked_R = L3,
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        
        return(xcos)
      })
      toc()
    }else{
      set.seed(100)
      tic("XCos")
      res_xcos<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        pval <- mat_pval[,i]
        names(fc) <- gene_names
        S <-get_deg(fc, abs_thresh = 1.5, pval = pval, use.pval = TRUE, thresh_pval = 0.05)
        xcos <- XCos_score(S = S, gene_names = gene_names,
                           pert_names = pert_names,list_ranked_R = L3,
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        return(xcos)
      })
      toc()
    }
  }else{
    set.seed(100)
    tic("Xcos")
    res_xcos<- purrr::map(seq(1,ncol(mat_queries)),function(i){
      fc <- mat_queries[,i]
      names(fc) <- gene_names
      S <-get_deg(fc, use.top = TRUE, topN= topN)
      xcos <- XCos_score(S = S, gene_names = gene_names,
                         pert_names = pert_names,list_ranked_R = L3,
                         ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
      
      return(xcos)
    })
    toc()
  }
}else{
  if(is.na(topN)){
    if(is.null(mat_pval)){
      set.seed(100)
      tic("Xcos")
      res_xcos<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        names(fc) <- gene_names
        sel <- which(params$ind_base==i)
        S <-get_deg(fc, abs_thresh = 1.5, pval = NA, use.pval = FALSE, thresh_pval = 0.05)
        xcos <- XCos_score(S = S, gene_names = gene_names,
                           pert_names = pert_names[sel],list_ranked_R = L3[sel],
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        
        return(xcos)
      })
      toc()
    }else{
      set.seed(100)
      tic("XCos")
      res_xcos<- purrr::map(seq(1,ncol(mat_queries)),function(i){
        fc <- mat_queries[,i]
        pval <- mat_pval[,i]
        names(fc) <- gene_names
        sel <- which(params$ind_base==i)
        S <-get_deg(fc, abs_thresh = 1.5, pval = pval, use.pval = TRUE, thresh_pval = 0.05)
        xcos <- XCos_score(S = S, gene_names = gene_names,
                           pert_names = pert_names[sel],list_ranked_R = L3[sel],
                           ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
        return(xcos)
      })
      toc()
    }
  }else{
    set.seed(100)
    tic("Xcos")
    res_xcos<- purrr::map(seq(1,ncol(mat_queries)),function(i){
      fc <- mat_queries[,i]
      names(fc) <- gene_names
      sel <- which(params$ind_base==i)
      S <-get_deg(fc, use.top=TRUE, topN=topN)
      xcos <- XCos_score(S = S, gene_names = gene_names,
                         pert_names = pert_names[sel],list_ranked_R = L3[sel],
                         ncpus = ncpus,permuted_pval = permuted_pval,permute_nb = 100)
      
      return(xcos)
    })
    toc()
  }
}

saveRDS(res_xcos,paste0(save_path,"_",topN,".rds"))
