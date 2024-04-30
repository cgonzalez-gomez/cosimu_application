########### SIGNATURES SIMULATION : REAL DATA INFERED PARAMETRIZATION ##########

# Imports
library(purrr)
source("./paper_scripts/simulation_pipeline.R")

path <- "./paper_scripts/RDS_symmetrical" 
save_param <- paste0(path,"df_param.RDS")
base_seed <- 19

# Load distribution parameters inferred from real data
dist <- cosimu::load_dist("./data/real_symmetrical_dist.yml")

# Random number generator
RNGkind("L'Ecuyer-CMRG")


# Parametrization for the primary and secondary signatures
## Primary
nr_noise <- c(0.02,0.1)
list_pDEG <- list(c(p_up=0.1, p_down=0.1),
                  c(p_up=0.05, p_down=0.15))

p_nb_mol_ctrl <- readRDS("data/ptt_real.RDS")
p_nb_mol_ctrl <- p_nb_mol_ctrl[!is.na(p_nb_mol_ctrl)]
gene_names <- names(p_nb_mol_ctrl)
nb_genes = length(gene_names)

ncpus <- 30
nb_tech_rep <- 20

ind_param <- map(rep(list_pDEG,each=nb_tech_rep),function(pDEG){
  list(p_up= pDEG[1],p_down=pDEG[2],prop_sm_up=dist$prop_sm_up,
       prop_sm_down=dist$prop_sm_down,
       qf_vect_up=dist$qf_vect_up,qf_vect_down=dist$qf_vect_down,
       qf_vect_nr=dist$qf_vect_nr,up_means= dist$up_means,
       nr_means=dist$nr_means,down_means= dist$down_means)
})

## Secondary
submod_transition <- "cop"
copula_submod <- "Frank"
rho_submod <- 0.9
eps_submod <- 1e-3
optim_method_submod <- "Brent"
proba_transition <- "deterministic"

set.seed(100)
cs <- c(rbeta(300,5,1.2),rbeta(200,1.5,3))

dep_param_df <- expand.grid(cs,nr_noise)
colnames(dep_param_df) <- c("cs","nr_noise")

dep_param = map(seq(1,nrow(dep_param_df)),function(i){ 
  list(mod_transition="default",
       submod_transition=submod_transition,
       qf_vect_up=dist$qf_vect_up,
       qf_vect_down=dist$qf_vect_down,
       qf_vect_nr=dist$qf_vect_nr,
       expected_cs= dep_param_df$cs[i],
       transition_mat=NULL,nr_noise= dep_param_df$nr_noise[i],
       prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
       proba_transition=proba_transition, 
       rho_submod= rho_submod, copula_submod = copula_submod,
       eps_submod = eps_submod,
       optim_method_submod = optim_method_submod,up_means= dist$up_means,
       nr_means=dist$nr_means,down_means= dist$down_means)
}
)

idx_ind_base <- rep(seq(1,length(ind_param)), each=length(dep_param))
idx_dep_base <- rep(seq(1,length(dep_param)),times=length(ind_param))

## Global parametrization data frame
gb_param_df <- data.frame(ind_base= idx_ind_base,
                          dep_base=idx_dep_base)
gb_param_df <- cbind(gb_param_df,dep_param_df, 
                     do.call(
                       rbind,rep(list_pDEG, each=nb_tech_rep*length(dep_param))
                     )
)
gb_param_df$dep_id <- paste0("dep_",seq(1,nrow(gb_param_df)))
gb_param_df$ind_id <- paste0("ind_",gb_param_df$ind_base)
saveRDS(gb_param_df,save_param)

# Parametrization for the replicates interconnected signatures
cs_inter <- 0.9
nr_noise_inter <- 0.02
submod_transition_inter <- "cop"
copula_submod_inter <- "Frank"
rho_submod_inter <- 0.9
eps_submod_inter <- 1e-3
optim_method_submod_inter <- "Brent"
proba_transition_inter <- "cop"
copula_prob_inter <- "deterministic"
inter_param = list(mod_transition="default",
                   submod_transition=submod_transition_inter,
                   qf_vect_up=dist$qf_vect_up,
                   qf_vect_down=dist$qf_vect_down,
                   qf_vect_nr=dist$qf_vect_nr,
                   expected_cs= cs_inter,
                   transition_mat=NA,nr_noise= nr_noise_inter,
                   prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                   proba_transition=proba_transition_inter, 
                   rho_submod= rho_submod_inter, copula_submod = copula_submod_inter,
                   eps_submod = eps_submod_inter,
                   optim_method_submod = optim_method_submod_inter)

# Seeds 
set.seed(base_seed)
list_seeds <- map(unique(gb_param_df$ind_base),~sample(seq(1,5000),7,replace = F))

for(i in unique(gb_param_df$ind_base)){
  dir_name <- paste0(path,"Q_",i)
  dir.create(dir_name, showWarnings = FALSE)
  suppressWarnings(simulation_pipeline(ncpus=ncpus, nb_genes=nb_genes, 
                                       ind_param=ind_param[i],
                                       gene_names=gene_names,
                                       dep_param = dep_param,
                                       p_nb_mol=p_nb_mol_ctrl,
                                       inter_param = inter_param,
                                       seeds=list_seeds[[i]],
                                       save_lfc_cosimu = 
                                         paste0(dir_name,"/cosimu.RDS"),
                                       save_DEA = 
                                         paste0(dir_name,"/deseq2.RDS"),
                                       ncpus_secondary = 3))
  
}

