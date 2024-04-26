# Parametrization to generate 100 interconnected signatures with the distribution
# parameters inferred from real log2 fold-change vector obtain by analyzing GSE185985 data.
# Counts were analyzed with the default DESeq2 parametrization, the expression vector 
# corresponds to the comparison between DMSO and TLN468 treatment.

# Imports
#library(cosimu)
library(purrr)
library(yaml)


# Random number generator
RNGkind("L'Ecuyer-CMRG")


# Generate YAML file with the dist parametrization (only needs to be run ones)
# lfc_vec <- read.csv("./data/real_lfc_data.csv",header = T)$log2FoldChange
# cosimu::infer_dist_param(lfc_vec,save_yaml = "./data/real_asymmetrical_dist.yml")

# Load the distribution parametrization (not need to load the functions, so it's
# not necessary to call load_dist_param)
dist <- read_yaml("./data/real_asymmetrical_dist.yml")

# Primary signature parametrization 
p_nb_mol_ctrl <- readRDS("./data/ptt_real.RDS")
p_nb_mol_ctrl <- p_nb_mol_ctrl[!is.na(p_nb_mol_ctrl)]
gene_names <- names(p_nb_mol_ctrl)
base_expression <- TRUE
nb_tech_rep <- 1

## Primary
nb_ent <- length(gene_names)
pDEG <- 0.2
list_pDEG <- list(c(p_up=pDEG*dist$alpha, p_down=pDEG*(1-dist$alpha)))

## Secondary
nr_noise <- 0.02
mod_transition <- "default"
transition_mat <- NULL
submod_transition <- "cop"
copula_submod <- "Frank"
rho_submod <- 0.9
eps_submod <- 1e-3
optim_method_submod <- "Brent"
proba_transition <- "deterministic"
copula_prob <- NULL
theta_prob <- NULL
nbins_prob <- NULL

set.seed(100)
connectivity_score <- c(rbeta(80,5,1.2),rbeta(20,1.5,3))

# YAML
write_yaml(list("dist"=dist,
                "signatures"=list(
                  "initial_base"=p_nb_mol_ctrl,
                  "nb_tech_rep"=nb_tech_rep,
                  "base_expression"=base_expression,
                  "entity_id"=gene_names,
                  "primary"=list(
                    "nb_ent" =nb_ent,
                    "list_pDEG" = list_pDEG
                  ),
                  "secondary"=list(
                    "mod_transition"=mod_transition,
                    "transition_mat"=transition_mat,
                    "nr_noise"=nr_noise,
                    "connectivity_score"=connectivity_score,
                    "submod_transition"=submod_transition,
                    "copula_submod"=copula_submod,
                    "rho_submod"=rho_submod,
                    "eps_submod"=eps_submod,
                    "optim_method_submod"=optim_method_submod,
                    "proba_transition"=proba_transition,
                    "copula_prob"=copula_prob,
                    "theta_prob"=theta_prob,
                    "nbins_prob"=nbins_prob
                  )
                )
              ), file = "./data/example_param.yaml")
