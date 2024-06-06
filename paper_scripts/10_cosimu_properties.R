library(parallel)
library(cosimu)
library(purrr)
library(tictoc)

############################ SIMULATION ############################ 

scale_r = 0.3
relative_prop = c(0.4,0.35,0.25,0.1) 
shapes <- c(4,7,12,15)

dist <- list(prop_sm_up = relative_prop/sum(relative_prop),
             prop_sm_down = relative_prop/sum(relative_prop),
             qf_vect_up = purrr::map(
               shapes,
               function(shape) {
                 local(function(x) {
                   p <- parent.env(environment())
                   p$shape <- shape
                   p$scale_r <- scale_r
                   qgamma(x,shape,scale = scale_r)
                 })
               }
             ),
             qf_vect_nr = list(
               function(x) qnorm(x,0,0.15)
             ),
             qf_vect_down = purrr::map(
               shapes, 
               function(shape) {
                 local(function(x) {
                   p <- parent.env(environment())
                   p$shape <- shape
                   p$scale_r <- scale_r
                   -qgamma(x,shape,scale = scale_r)
                 })
               }
             ))

cs <- c(0,0.5,0.8,1)
score <- "pearson" 
nr_noise <- c(0.02,0.05,0.1,0.15)
pDEG <- 0.1
alpha <- 0.5

submod_transition <- c("independent","stochastic","cop","deterministic")
copula_submod <- "Frank"
rho_submod <- 0.9
eps_submod <- 1e-3
optim_method_submod <- "Brent"

proba_transition <- c("independent","cop","deterministic")
copula_prob <- "Frank"
theta_prob <- 10
nbins_prob <- 100 

ref_noise <- 0.02
ref_sm <- "deterministic"
ref_proba <- "deterministic"

param_df <- expand.grid(cs,nr_noise,ref_sm,ref_proba)
param_df <- rbind(param_df,expand.grid(cs,ref_noise,submod_transition,ref_proba))
param_df <- rbind(param_df,expand.grid(cs,ref_noise,ref_sm,proba_transition))

colnames(param_df) <- c("cs","nr_noise","submod_transition","proba_transition")
nb_ent <- 10000
gene_names <- NULL

nb_tech_rep <- 100
tech_rep <- rep(seq(1,nrow(param_df)),each=nb_tech_rep)
param_df<- param_df[tech_rep,]
param_df$tech_rep <- tech_rep


estimation <- function(nb_ent,dist,submod_transition,pDEG,
                       proba_transition,expected_cs,nr_noise,alpha,score){
  p_up <- pDEG*alpha
  p_down <- pDEG*(1-alpha)
  
  prim_sig <- PrimarySignatureObj$new(nb_ent = nb_ent,p_up = p_up,
                                      p_down = p_down,
                                      prop_sm_up = dist$prop_sm_up,
                                      prop_sm_down = dist$prop_sm_down,
                                      qf_vect_up = dist$qf_vect_up,
                                      qf_vect_nr = dist$qf_vect_nr,
                                      qf_vect_down = dist$qf_vect_down,
                                      entity_id = gene_names,
                                      base_expression = FALSE)
  sec_sig <- SecondarySignatureObj$new(
    prim_sig_obj = prim_sig, mod_transition = "default",
    submod_transition = submod_transition,
    proba_transition = proba_transition,
    qf_vect_up = dist$qf_vect_up, qf_vect_nr = dist$qf_vect_nr,
    qf_vect_down = dist$qf_vect_down,
    connectivity_score = expected_cs,
    nr_noise = nr_noise, 
    prop_sm_up = dist$prop_sm_up,
    prop_sm_down = dist$prop_sm_down,
    rho_submod = rho_submod,
    copula_submod = copula_submod,
    eps_submod = eps_submod,
    optim_method_submod = optim_method_submod,
    nbins_prob = nbins_prob,
    theta_prob = theta_prob,
    copula_prob = copula_prob,
    base_expression = FALSE)
  return(cor(prim_sig$get_lfc_vect(),sec_sig$get_lfc_vect(),method=as.character(score)))
}

nb_cores <- 7
set.seed(123)
tictoc::tic()
res<-mclapply(seq(1,nrow(param_df)), function(i){ 
  estimation(nb_ent = nb_ent,dist = dist,
             submod_transition = param_df$submod_transition[i],
             proba_transition = param_df$proba_transition[i],
             pDEG = pDEG,alpha=alpha,
             expected_cs = param_df$cs[i],nr_noise = param_df$nr_noise[i],
             score = score)
},mc.cores = nb_cores)
tictoc::toc()

# saveRDS(res, path)
# saveRDS(param_df, path)
############################ PLOT ############################ 
library(ggplot2)
data_final <- res #readRDS(path)
param_final <- param_df #readRDS(path)

colnames(param_final) <- c("cs","nr_noise","submod_transition","proba_transition","tech_rep")
df <- cbind(param_final,cor=unlist(data_final))
df$nr_noise <- as.factor(df$nr_noise)
df$submod_transition <- factor(df$submod_transition,levels = c("deterministic","cop","stochastic","independent"))
df$proba_transition <- factor(df$proba_transition,levels = c("deterministic","cop","independent"))

ref_noise <- 0.02
ref_sm <- "deterministic"
ref_proba <- "deterministic"

s = 25 #text size
### modality
sel_mod<-which(df$submod_transition==ref_sm  & df$proba_transition==ref_proba)

m <- ggplot(df[sel_mod,],aes(x=nr_noise,y=cor,color=nr_noise))+
  geom_boxplot()+
  ggh4x::facet_nested("connectivity score" + factor(cs,levels=c(1,0.8,0.5,0)) ~ .,
                      scales = "free")+
  theme(axis.text.x = element_text(size = s),
        axis.text.y = element_text(size = s),
        axis.title = element_text(size = s),
        strip.text = element_text(size = s),
        legend.position = "none")+
  scale_y_continuous(n.breaks = 4)+
  xlab(expression("Modality transition noise factor" ~ (gamma)))+
  ylab("Pearson's correlation")+
  scale_color_manual(values=c("#B50003","#F46D43","#74ADD1","#244B7E"))
m

### submodality
sel_submod<-which(df$nr_noise==ref_noise 
                  & df$proba_transition==ref_proba)

sm <- ggplot(df[sel_submod,],aes(x=submod_transition,y=cor,color=submod_transition))+
  geom_boxplot()+
  ggh4x::facet_nested("connectivity score" + factor(cs,levels=c(1,0.8,0.5,0)) ~ .,
                      scales = "free")+
  theme(axis.text.x = element_text(size = s),
        axis.text.y = element_text(size = s),
        axis.title = element_text(size = s),
        strip.text = element_text(size = s),
        legend.position = "none")+
  scale_y_continuous(n.breaks = 4)+
  xlab("Sub-modality transition")+
  ylab("Pearson's correlation")+
  scale_x_discrete(labels=c("det.","cop.","stoch.","ind."))+
  scale_color_manual(values=c("#B50003","#F46D43","#74ADD1","#244B7E"))
sm

### probability 
sel_proba<-which(df$nr_noise==ref_noise 
                 & df$submod_transition==ref_sm)

p <- ggplot(df[sel_proba,],aes(x=proba_transition,y=cor,color=proba_transition))+
  geom_boxplot()+
  ggh4x::facet_nested("connectivity score" + factor(cs,levels=c(1,0.8,0.5,0)) ~ .,
                      scales = "free")+
  theme(axis.text.x = element_text(size = s),
        axis.text.y = element_text(size = s),
        axis.title = element_text(size = s),
        strip.text = element_text(size = s),
        legend.position = "none")+
  scale_y_continuous(n.breaks = 4)+
  xlab("Probability transition")+
  ylab("Pearson's correlation")+
  scale_x_discrete(labels=c("det.","cop.","ind."))+
  scale_color_manual(values=c("#B50003","#74ADD1","#244B7E"))

p

final_plt <- ggpubr::ggarrange(m,sm,p, ncol = 3,nrow = 1)
final_plt

