#' Generation of the read counts from the fold-change signatures. This function 
#' was inspired from the first part of the function `simulate_experiment` from the
#' package newpolyester. 
#' (Frazee AC, Jaffe AE, Langmead B, Leek JT. 
#' Polyester: simulating RNA-seq datasets with differential transcript expression. 
#' Bioinformatics. 2015 Sep 1;31(17):2778-84. doi: 10.1093/bioinformatics/btv272.
#' Epub 2015 Apr 28. PMID: 25926345; PMCID: PMC4635655.) 
#' @param num_reps: number of replicates for each condition
#' @param p_nb_mol: proportional to true number of molecules associated to each
#' gene in its a basal expression level (untreated) (inspired from DESeq2 
#' proportional to true number of fragments, default 300) 
#' @param size: optional size for the NB model
#' @param size_factor: optional factor that set the size of the NB to mean/size_factor
#' @param fold_changes: fold-changes values
#' @param num_reads: expected number of reads in each experiment
#' @param TP_analysis: boolean set to TRUE if a 3'-RNASeq analysis parametrization
#' is considered. This impacts the influence of the gene or transcripts length 
#' in the analysis, that is not take into account in a 3' analysis.
#' @param length_ent: gene or transcripts lengths.
#' @param nb_ent: number of genes or transcripts.
#' @param names_ent: names of the genes or transcripts.
#' @param seed: random seed.
read_counts_gen <- function(num_reps=c(3,3), p_nb_mol=300,size=NULL,
                            size_factor = 3,fold_changes,num_reads=1000000,
                            TP_analysis=FALSE,length_ent = NULL,
                            nb_ent = NULL, names_ent = NULL, seed = NULL){
  # Seeds
  if(!is.null(seed)){
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }
  ###########
  # Validation 
  if(length(num_reps) == 1 & !is.null(dim(fold_changes))){
    num_reps = rep(num_reps, ncol(fold_changes))
  }
  
  stopifnot(is.matrix(fold_changes))
  if(ncol(fold_changes) != length(num_reps)){
    stop(.makepretty('wrong number of columns in fold change matrix:
                need same number of columns as number of groups.'))
  }
  if(nrow(fold_changes) != nb_ent){
    stop(.makepretty('wrong number of rows in fold change matrix: need
                same number of rows as number of simulated entities.'))
  }

  ##########
  # Basemeans 
  ptt = p_nb_mol * fold_changes
  if(TP_analysis) {
    Lg = 1
  }else{
    Lg=  length_ent
  }
  S=colSums(matrix(ptt*Lg,ncol=length(num_reps)))
  basemeans =  sweep(x=as.matrix(ptt*Lg,ncol=length(num_reps)),MARGIN = 2,STATS =S,FUN = "/")*num_reads
  ##########
  # Size factor 
  if(is.null(size)){
    size = basemeans / size_factor
    size[size==0]=0.0001
  }else if(class(size) == 'numeric'){
    if(any(size==0)){
      size[size==0]=0.0001
    }
    size = matrix(size, nrow=nrow(basemeans), ncol=ncol(basemeans))
  }else if(inherits(size, 'matrix')){
    stopifnot(nrow(size) == nrow(basemeans))
    stopifnot(ncol(size) == ncol(basemeans))
  }else{
    stop('size must be a number, numeric vector, or matrix.')
  }
  ###########
  # Readmat
  group_ids = rep(1:length(num_reps), times=num_reps)
  numreadsList = vector("list", sum(num_reps))
  
  if(inherits(size, 'matrix')){
    numreadsList = lapply(1:sum(num_reps), function(i){
      group_id = group_ids[i]
      return(rnbinom(n = nb_ent, mu = basemeans[,group_id], size = size[,group_id]))
    })
  }else if(class(size) == 'numeric'){
    numreadsList = lapply(1:sum(num_reps), function(i){
      group_id = group_ids[i]
      return(rnbinom(n = nb_ent, mu = basemeans[,group_id], size = size[group_id]))
    })
  }
  readmat = matrix(unlist(numreadsList), ncol=sum(num_reps))
  ###########
  # Names
  if(is.null(colnames(fold_changes))){
    colnames(readmat) <- paste0("sample_", seq(1, sum(num_reps)))
  }else{
    idx <- unlist(purrr::map(num_reps, ~seq_len(.x)))
    colnames(readmat) <- paste0(rep(colnames(fold_changes), num_reps), "_",idx)
  }
  rownames(readmat) <- names_ent
  ###########  
  return(readmat)
}


