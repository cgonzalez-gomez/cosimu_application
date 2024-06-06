# Simulation of 100 interconnected signatures based on a YAML parametrization,
# without biological replicates
library(cosimu)
# Load the parametrization usinng the function load_simple_param from cosimu
parameters <- load_simple_param("./data/example_param.yaml") # 
entity_id <- parameters$entity_id
base_expression <- parameters$base_expression
initial_base <- parameters$initial_base
nb_ent <- length(initial_base)
primary_param <- parameters$primary_param
secondary_param <- parameters$secondary_param
idx_p <- parameters$gb_param_df$primary_base
idx_s <- parameters$gb_param_df$secondary_base

ncpus <- 1
set.seed(100)
parallel::mc.reset.stream()
list_primary_sign <- parallel::mclapply(primary_param, function(param){
  res <- PrimarySignatureObj$new(nb_ent = nb_ent,p_up = param$p_up,
                                 p_down = param$p_down,
                                 prop_sm_up = param$prop_sm_up,
                                 prop_sm_down = param$prop_sm_down,
                                 qf_vect_up = param$qf_vect_up,
                                 qf_vect_nr = param$qf_vect_nr,
                                 qf_vect_down = param$qf_vect_down,
                                 entity_id = entity_id,
                                 base_expression = base_expression,
                                 initial_base = initial_base,
                                 up_means = param$up_means,
                                 nr_means = param$nr_means,
                                 down_means = param$down_means)
  return(res)
}, mc.cores = min(length(primary_param), ncpus), mc.set.seed = TRUE)

# Treatment 2
set.seed(200)
parallel::mc.reset.stream()
list_secondary_sign <- parallel::mclapply(
  seq(1,(length(secondary_param)*length(list_primary_sign))),function(i){
    param <- secondary_param[[idx_s[i]]]
    prim_sig <- list_primary_sign[[idx_p[i]]]
    obj <- SecondarySignatureObj$new(
      prim_sig_obj = prim_sig, mod_transition = param$mod_transition,
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
      base_expression = base_expression,
      up_means = param$up_means,nr_means = param$nr_means,
      down_means = param$down_means,
      ncpus = ncpus)
    return(obj)
  }, mc.cores = min(length(secondary_param), ncpus), mc.set.seed = TRUE
)

# Save the signatures matrix 
df_lfc <- cbind(purrr::map_dfc(list_primary_sign, ~ .x$get_lfc_vect()),
                    purrr::map_dfc(list_secondary_sign, ~.x$get_lfc_vect()))
colnames(df_lfc) <- c(paste0("primary_",seq(1,length(list_primary_sign))),
                          paste0("secondary_",seq(1,length(list_secondary_sign))))

rownames(df_lfc) <- entity_id