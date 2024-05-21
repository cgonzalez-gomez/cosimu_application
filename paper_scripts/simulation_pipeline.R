# PIPELINE COSIMU - READ COUNTS - DESEQ2
library(cosimu)
library(purrr)
library(DESeq2)
library(parallel)
library(tictoc)

# Inputs
#' @param ncpus: Number of CPUs to be used in parallel.
#' @param ncpus_secondary: The maximum number of CPUs to be used in parallel during
#' the secondary signatures simulation is ncpus_secondary x ncpus.
#' @param seeds: vector of random seeds to be used.
#' @param nb_genes: Length of the gene signature.
#' @param gene_names: Character vector with the genes names.
#' @param ind_param: List of named lists for the parametrization of the Primary
#' signatures. Each list contains the parameter of the `PrimarySignatureObj` of
#' cosimu.
#' @param inter_param : parametrization for the interaction between experimental
#' replicates named vector with the parameters of the `SecondarySignatureObj` of
#' cosimu.
#' @param dep_param (optional): List of named lists for the parametrization of the 
#' Secondary signatures (default NA, generates only primary vectors). Each list 
#' contains the parameters of the `SecondarySignatureObj` of cosimu.
#' @param save_lfc_cosimu: Path to save the log-fold change matrix.
#' @param add_low: boolean set to TRUE to include low values in the basal
#' signature without generating an error. This values will be included in the
#' read counts generation process but not in the bias correction.
#' @param basemeans_thr: threshold to filter low values from the basal signature 
#' @param p_nb_mol : proportional to true number of molecules (ppt) associated 
#' to each gene in its a basal expression level (untreated) (inspired from
#' DESeq2 proportional to true number of fragments).
#' @param num_reps: (int) number of expected experimental replicates.
#' @param num_reads: (int) Expected number of reads in each experiment 
#' (inherited from `read_count_gen.R`).
#' @param TP_analysis: Boolean set to TRUE if a 3'-RNASeq analysis
#' parametrization is considered (inherited from `cosimu::read_count_gen`)
#' @param size_factor: (inherited from `cosimu::read_count_gen`)
#' @param length_ent : gene or transcripts lengths (inherited from 
#' `read_count_gen.R`)
#' @param save_counts: Path to save the simulated read counts matrix.
#' @param DE_analysis: Boolean set to TRUE to run DESeq2 analysis.
#' @param return_res: Boolean set to TRUE to return a list with the following
#' elements: 
#'  - signature objects: primary, secondary signatures from
#'  both treatments
#'  - corrected ppt (after inference bias correction) for both treatments
#'  - list of ppt after the introduction of the biological variability
#'  - count matrix
#'  - DEA results
simulation_pipeline <- function(ncpus, nb_genes, ind_param,inter_param,
  gene_names,p_nb_mol,dep_param=NULL,
  seeds = c(123,147,78,300,400,596), num_reps = 3, save_lfc_cosimu = "./lfc_cosimu.RDS", 
  num_reads=4300000, TP_analysis=TRUE, save_counts = NULL,
  DE_analysis = TRUE ,save_DEA = "./DEA.RDS", return_res = FALSE,
  correct_BE = TRUE,sigma_bio=0.5,add_low=TRUE,basemeans_thr = 1,
  size_factor = 3, size = NULL, length_ent =NULL, paired = TRUE,
  ncpus_secondary = 1){
  
  if(add_low){
    full_p_nb_mol <- p_nb_mol
    full_nb_genes <- nb_genes
    full_names <- gene_names
    subset <- which(p_nb_mol > basemeans_thr & !is.na(p_nb_mol))
    p_nb_mol <- p_nb_mol[subset]
    nb_genes <- length(subset)
    gene_names <- gene_names[subset]
  }
  
  RNGkind("L'Ecuyer-CMRG")
  
  if(length(seeds)<6){
    seeds <- rep(seeds,length.out = 6)
  }

  if(return_res){
    glob_res <- list()
  }
  p_nb_mol <- p_nb_mol[gene_names]
  ################################### COSIMU ###################################
  set.seed(seeds[1])
  parallel::mc.reset.stream()
  list_ind_base <- parallel::mclapply(ind_param, function(param){
    res <- PrimarySignatureObj$new(nb_ent = nb_genes,p_up = param$p_up,
                                   p_down = param$p_down,
                                   prop_sm_up = param$prop_sm_up,
                                   prop_sm_down = param$prop_sm_down,
                                   qf_vect_up = param$qf_vect_up,
                                   qf_vect_nr = param$qf_vect_nr,
                                   qf_vect_down = param$qf_vect_down,
                                   entity_id = gene_names,
                                   base_expression = correct_BE,
                                   initial_base = p_nb_mol,
                                   up_means = param$up_means,
                                   nr_means = param$nr_means,
                                   down_means = param$down_means)
    return(res)
  }, mc.cores = min(length(ind_param), ncpus), mc.set.seed = TRUE)
  
  # Treatment 1
  tic("Signatures treatment #1 : secondary")
  list_ind <- parallel::mclapply(rep(seq(1,length(list_ind_base)),
                                     each=num_reps), function(i){
    base <- list_ind_base[[i]]
    param <- ind_param[[i]]
    obj <- SecondarySignatureObj$new(
      prim_sig_obj = base, mod_transition = inter_param$mod_transition,
      submod_transition = inter_param$submod_transition,
      proba_transition = inter_param$proba_transition,
      qf_vect_up = inter_param$qf_vect_up, qf_vect_nr = inter_param$qf_vect_nr,
      qf_vect_down = inter_param$qf_vect_down,
      connectivity_score = inter_param$connectivity_score,
      transition_mat = inter_param$transition_mat,
      nr_noise = inter_param$nr_noise, prop_sm_up = inter_param$prop_sm_up,
      prop_sm_down = inter_param$prop_sm_down,
      rho_submod = inter_param$rho_submod,
      copula_submod = inter_param$copula_submod,
      eps_submod = inter_param$eps_submod,
      optim_method_submod = inter_param$optim_method_submod,
      nbins_prob = inter_param$nbins_prob,
      theta_prob = inter_param$theta_prob,
      copula_prob = inter_param$copula_prob, 
      base_expression = correct_BE,
      up_means = param$up_means,nr_means = param$nr_means,
      down_means = param$down_means,
      ncpus = ncpus)
    return(obj)
  }, mc.cores = min((length(list_ind_base)*num_reps), ncpus_secondary),
  mc.set.seed = TRUE )
  toc()
  
  # Treatment 2
  if(!is.null(dep_param)){
    set.seed(seeds[2])
    parallel::mc.reset.stream()
    idx_ind_base <- rep(seq(1,length(list_ind_base)), each=length(dep_param))
    idx_param <- rep(seq(1,length(dep_param)),times=length(list_ind_base))
    tic("Signatures Treatment #2 : primary")
    list_dep_base <- parallel::mclapply(
      seq(1,(length(dep_param)*length(list_ind_base))),function(i){
        param <- dep_param[[idx_param[i]]]
        ind_base <- list_ind_base[[idx_ind_base[i]]]
        obj <- SecondarySignatureObj$new(
          prim_sig_obj = ind_base, mod_transition = param$mod_transition,
          submod_transition = param$submod_transition,
          proba_transition = param$proba_transition,
          qf_vect_up = param$qf_vect_up, qf_vect_nr = param$qf_vect_nr,
          qf_vect_down = param$qf_vect_down,
          connectivity_score = param$connectivity_score,
          transition_mat = param$transition_mat,
          nr_noise = param$nr_noise, prop_sm_up = param$prop_sm_up,
          prop_sm_down = param$prop_sm_down,
          rho_submod = param$rho_submod, copula_submod = param$copula_submod,
          eps_submod = param$eps_submod,
          optim_method_submod = param$optim_method_submod,
          nbins_prob = param$nbins_prob,
          theta_prob = param$theta_prob,
          copula_prob = param$copula_prob,
          base_expression = correct_BE,
          up_means = param$up_means,nr_means = param$nr_means,
          down_means = param$down_means,
          ncpus = ncpus)
        return(obj)
      }, mc.cores = min(length(dep_param), ncpus_secondary), mc.set.seed = TRUE
    )
    toc()
    tic("Signatures Treatment #2 : secondary")
    list_dep <- parallel::mclapply(
      rep(list_dep_base, each=num_reps), function(base){
        obj <- SecondarySignatureObj$new(
          prim_sig_obj = base, mod_transition = inter_param$mod_transition,
          submod_transition = inter_param$submod_transition,
          proba_transition = inter_param$proba_transition,
          qf_vect_up = inter_param$qf_vect_up,
          qf_vect_nr = inter_param$qf_vect_nr,
          qf_vect_down = inter_param$qf_vect_down,
          connectivity_score = inter_param$connectivity_score,
          transition_mat = inter_param$transition_mat,
          nr_noise = inter_param$nr_noise, prop_sm_up = inter_param$prop_sm_up,
          prop_sm_down = inter_param$prop_sm_down,
          rho_submod = inter_param$rho_submod,
          copula_submod = inter_param$copula_submod,
          eps_submod = inter_param$eps_submod,
          optim_method_submod = inter_param$optim_method_submod,
          nbins_prob = inter_param$nbins_prob,
          theta_prob = inter_param$theta_prob,
          copula_prob = inter_param$copula_prob, ncpus = ncpus)
        return(list(lfc=obj$get_lfc_vect()))
      }, mc.cores = min((length(list_ind_base)*num_reps), ncpus_secondary),
      mc.set.seed = TRUE
    )
    toc()
  }
  
  if(!is.null(save_lfc_cosimu)) { # save signatures matrix
    if(!is.null(dep_param)){
      df_fc_gene <- cbind(purrr::map_dfc(list_ind, ~ .x$get_lfc_vect()),
                          purrr::map_dfc(list_dep, ~.x$lfc))
      colnames(df_fc_gene) <- c(paste0("ind_",seq(1,length(list_ind))),
                                paste0("dep_",seq(1,length(list_dep))))
    }else{
      df_fc_gene <- purrr::map_dfc(list_ind, ~ .x$get_lfc_vect())
      colnames(df_fc_gene) <- paste0("ind_",seq(1,length(list_ind)))
    }
    rownames(df_fc_gene) <- gene_names
    saveRDS(df_fc_gene, save_lfc_cosimu)
    df_fc_gene <- 2^df_fc_gene # transform log2 fold-change into fold-change
  }else{
    # transform log2 fold-change into fold-change
    if(!is.null(dep_param)){
      df_fc_gene <- cbind(purrr::map_dfc(list_ind, ~ 2^.x$get_lfc_vect()),
                          purrr::map_dfc(list_dep, ~2^.x$lfc))
      colnames(df_fc_gene) <- c(paste0("ind_",seq(1,length(list_ind))),
                                paste0("dep_",seq(1,length(list_dep))))
    }else{
      df_fc_gene <- purrr::map_dfc(list_ind, ~ 2^.x$get_lfc_vect())
      colnames(df_fc_gene) <- paste0("ind_",seq(1,length(list_ind)))
    }
    rownames(df_fc_gene) <- gene_names
  }
  if(return_res){
    # Primary signatures from the treatments 1 & 2 
    glob_res$list_ind_base <- list_ind_base
    if(!is.null(dep_param)) glob_res$list_dep_base <- list_dep_base

    # Secondary signatures from the treatments 1 & 2
    glob_res$cosimu_fc <- df_fc_gene
    warning("The returned data associated to cosimu corresponds to the FC,
            while the saved data is the LFC")
  }
  ################## PROPORTIONAL TO TRUE (PTT) NUMBER OF MOL ##################
  ###################### FOR UNTREATED (control) SAMPLES #######################
  set.seed(seeds[3])
  parallel::mc.reset.stream()

  if(add_low){
    if(correct_BE){
        if(!is.null(dep_param)){
          bases <- map(c(list_ind_base,list_dep_base),
                       ~.x$get_base_expression())
        }
        else{
          bases <- map(list_ind_base, ~.x$get_base_expression())
        }
      
      # Corrected ppt values for treatments 1 & 2
      if(return_res) glob_res$bases <- bases 
      if(paired){ # biological variability (donnor)
        p_nb_mol_list <- 
          map(rep(seq(1,length(bases)),each=num_reps), function(i){
            tmp <- full_p_nb_mol
            tmp[subset] <- bases[[i]]*2^(rnorm(nb_genes,0,sigma_bio))
            return(tmp)
          })
      }else{
        p_nb_mol_list <- 
          map(rep(seq(1,length(bases)),each=num_reps),function(i){
            tmp <- full_p_nb_mol
            tmp[subset] <- bases[[i]]
            return(tmp)
          })
      }
     
    }else{
      if(paired){ # biological variability (donnor)
        p_nb_mol_list <- map(seq(1,num_reps), function(i){
          tmp <- full_p_nb_mol
          tmp[subset] <- p_nb_mol*2^(rnorm(nb_genes,0,sigma_bio))
          return(tmp)
        })
      }else{
        p_nb_mol_list <- map(seq(1,num_reps), function(i){
          tmp <- full_p_nb_mol
          tmp[subset] <- p_nb_mol
          return(tmp)
        })
      }
      
    }
    # Final ppt list (one per replicated)
    if(return_res) glob_res$p_nb_mol_list <- p_nb_mol_list
    
    low_fc <- data.frame(matrix(1,ncol=ncol(df_fc_gene),
                                nrow=length(full_p_nb_mol)-length(p_nb_mol)))
    rownames(low_fc) <- full_names[-subset]
    colnames(low_fc) <- colnames(df_fc_gene)
    df_fc_gene <- rbind(df_fc_gene, low_fc)
    df_fc_gene <- df_fc_gene[full_names,]
    
    # Reset full data 
    p_nb_mol <- full_p_nb_mol 
    nb_genes <-  full_nb_genes 
    gene_names <-  full_names 
    
  }else{
    if(correct_BE){

        if(!is.null(dep_param)){
          bases <- map(c(list_ind_base,list_dep_base),
                       ~.x$get_base_expression())
        }
        else{
          bases <- map(list_ind_base, ~.x$get_base_expression())
        }
      
      if(return_res) glob_res$bases <- bases
      if(paired){
        p_nb_mol_list <- 
          map(rep(seq(1,length(bases)),each=num_reps), function(i){
            bases[[i]]*2^(rnorm(nb_genes,0,sigma_bio))
          })
      }else{
        p_nb_mol_list <- 
          map(rep(seq(1,length(bases)),each=num_reps), function(i){
            bases[[i]]
          })
      }
      
    }else{
      if(paired){
        p_nb_mol_list <- map(seq(1,num_reps),
                             ~p_nb_mol*2^(rnorm(nb_genes,0,sigma_bio)))
      }else{
        p_nb_mol_list <- map(seq(1,num_reps), ~p_nb_mol)
      }
      
    }
    if(return_res) glob_res$p_nb_mol_list <- p_nb_mol_list
  }
  ########################### READ COUNTS GENERATION ########################### 
  tictoc::tic("Read Counts generation")
  set.seed(seeds[4])
  parallel::mc.reset.stream()
  # matrix of ones, not noise added for technical replicates
  df_base_fc <- matrix(1, nrow = nb_genes, ncol = num_reps,
                       dimnames = list(NULL,paste0("base_",seq_len(num_reps))))
  
  fc_mat <- as.matrix(df_fc_gene)
  c_names = colnames(fc_mat)
  d <- seq(1,ncol(fc_mat))
  chunks = split(d, ceiling(seq_along(d)/num_reps))
  names(chunks) <- NULL
  cts <- mclapply(chunks,function(chunk){
    df <-  
      purrr::map2_dfc(seq(1,num_reps),p_nb_mol_list[chunk], 
                      ~cosimu::read_counts_gen(p_nb_mol=.y,
                                       fold_changes = cbind(fc_mat[,chunk[.x]],
                                                            df_base_fc[,.x]),
                                       TP_analysis = TP_analysis,
                                       num_reps = c(1,1), num_reads = num_reads,
                                       nb_ent = nb_genes,names_ent =gene_names,
                                       size_factor = size_factor,
                                       length_ent = length_ent, size = size))
    
    df <- df[,c(seq(1,ncol(df),2),seq(2,ncol(df),2))]
    colnames(df) <- c(colnames(fc_mat[,chunk]),colnames(df_base_fc))
    return(as.matrix(df))
  },mc.cores = min(length(chunks),ncpus), mc.set.seed = TRUE )
  tictoc::toc()
  if(!is.null(save_counts)){
    saveRDS(cts, save_counts)   
  } 
  if(return_res) glob_res$cts <- cts 
  ###################### DIFFERENTIAL EXPRESSION ANALYSIS ######################
  tictoc::tic("DESeq2")
  if(DE_analysis){
    if (paired){
      design = "~ condition + donnor"
    }else{
      design = "~ condition"
    }
    deseq_res <- mclapply(1:length(cts), function(i){
      coldata <- data.frame(condition = rep(c("treated", "untreated"),
                                            each = num_reps),
                            donnor = paste0("D",seq(1,num_reps)))
      cts_df <- cts[[i]]
      rownames(coldata) <- colnames(cts_df)
      dds <- DESeqDataSetFromMatrix(countData = cts_df,
                                    colData = coldata,
                                    design = as.formula(design))
      dds$condition <- relevel(dds$condition, ref = "untreated")
      dds <- estimateSizeFactors(dds)
      dds <- DESeq(dds,quiet = TRUE)
      res <- results(dds, name = resultsNames(dds)[2],
                     independentFiltering = TRUE) 
      return(res)
    },mc.cores = min(length(cts),ncpus), mc.set.seed = TRUE )
    tictoc::toc()
    if(return_res){
      glob_res$deseq_DEA <- deseq_res
    }
    if(!is.null(save_DEA)){
      saveRDS(deseq_res, save_DEA)
    }
  }
  if(return_res) return(glob_res)
}


