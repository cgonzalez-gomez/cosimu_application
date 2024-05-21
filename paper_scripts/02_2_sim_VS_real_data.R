############# COMPARISON OF THE LOG2 FOLD CHANGE DISTRIBUTION ############# 
library(GEOquery)
library(stringr)
library(philentropy)
library(ggplot2)
library(purrr)
source('./paper_scripts/simulation_pipeline.R')
path_data <- "./downloaded_real_data/" #NEED UPDATE
path_data <- "~/Documents/signia_backup/real_data_for_corr/"
deseq_analysis<-function(cts,coldata,design){
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = as.formula(design))
  dds$condition <- relevel(dds$condition, ref = "untreated")
  norm_fact <- matrix(rep(1e+03, length(cts)), nrow = nrow(cts),
                      ncol = ncol(cts),dimnames = list(
                        row.names(cts),colnames(cts)))
  norm_fact <- norm_fact / exp(rowMeans(log(norm_fact)))
  dds <- estimateSizeFactors(dds, normMatrix = norm_fact)
  dds <- DESeq(dds,quiet = TRUE)
  res <- results(dds, name = resultsNames(dds)[2], independentFiltering = TRUE)
  return(res)
}

lim_padj <- 80

###### KL divergence between biological replicates ######
# Case 1 
data_GSE182024 <- getGEO(filename = paste0(path_data,"GSE182024_family.soft"))

sample_id <- data_GSE182024@header$sample_id
description <- map_chr(sample_id,~data_GSE182024@gsms[[.x]]@header$title)

cts_df <- read.table(file = paste0(path_data,'GSE182024_Raw_gene_counts_matrix.txt'), sep = '\t', header = TRUE)
cts_df <- cts_df[,description]

cts_df<-cts_df[rowSums(cts_df)!=0,]
characteristics <- map_dfr(description,function(x) {
  row <- str_split(x,"_")[[1]]
  ctrl <- "Control" %in% row
  TNFa <- "TNFa" %in% row
  return(c(ctrl=ctrl,IL4_IL13 = (!ctrl & !TNFa),TNFa=TNFa))
})

characteristics$description <- description

ctrl <- which(characteristics$ctrl==TRUE)
treat <- which(characteristics$IL4_IL13==TRUE)

coldata <- data.frame(condition = rep(c("untreated", "treated"), each = 3))
design = "~ condition"

res1 <- deseq_analysis(cts = cts_df[,c(ctrl[1:3],treat[1:3])],coldata = coldata, design = design)
res2 <- deseq_analysis(cts = cts_df[,c(ctrl[4:6],treat[4:6])],coldata = coldata, design = design)

lfc_rep1 <- res1$log2FoldChange
lfc_rep1 <-res1$log2FoldChange[-log10(res1$padj)<lim_padj]
lfc_rep1 <- lfc_rep1[!is.na(lfc_rep1)]

lfc_rep2 <- res2$log2FoldChange
lfc_rep2 <-res2$log2FoldChange[-log10(res2$padj)<lim_padj]
lfc_rep2 <- lfc_rep2[!is.na(lfc_rep2)]

h_real <- hist(lfc_rep1,breaks = seq(-7.1,17,0.01),freq = F)
h_real2 <- hist(lfc_rep2,breaks = seq(-7.1,17,0.01),freq = F)


kl1 <- KL(rbind(h_real2$counts/sum(h_real2$counts),h_real$counts/sum(h_real$counts)))
#kullback-leibler : 0.1264525

kl2 <- KL(rbind(h_real$counts/sum(h_real$counts),h_real2$counts/sum(h_real2$counts)))
# kullback-leibler : 0.1178246

mean_kl <- mean(c(kl1,kl2))
#mean_kl : 0.1221385


# Case 2
treat2 <- which(characteristics$IL4_IL13==FALSE & characteristics$TNFa==TRUE)

res2.1 <- deseq_analysis(cts = cts_df[,c(ctrl[1:3],treat2[1:3])],coldata = coldata, design = design)
res2.2 <- deseq_analysis(cts = cts_df[,c(ctrl[4:6],treat2[4:6])],coldata = coldata, design = design)

lfc_rep2.1 <- res2.1$log2FoldChange
lfc_rep2.1 <-res2.1$log2FoldChange[-log10(res2.1$padj)<lim_padj]
lfc_rep2.1 <- lfc_rep2.1[!is.na(lfc_rep2.1)]

lfc_rep2.2 <- res2.2$log2FoldChange
lfc_rep2.2 <-res2.2$log2FoldChange[-log10(res2.2$padj)<lim_padj]
lfc_rep2.2 <- lfc_rep2.2[!is.na(lfc_rep2.2)]

h_real2.1 <- hist(lfc_rep2.1,breaks = seq(-10.9,12.5,0.01),freq = F)
h_real2.2 <- hist(lfc_rep2.2,breaks = seq(-10.9,12.5,0.01),freq = F)


kl2.1 <- KL(rbind(h_real2.2$counts/sum(h_real2.2$counts),h_real2.1$counts/sum(h_real2.1$counts)))

kl2.2 <- KL(rbind(h_real2.1$counts/sum(h_real2.1$counts),h_real2.2$counts/sum(h_real2.2$counts)))

mean_kl2 <- mean(c(kl2.1,kl2.2))


###### KL divergence between simulated distribution and real distribution  ######
###### from which parameters where infered ######

# Load distribution parameters inferred from real data
dist <- cosimu::load_dist("./data/real_asymmetrical_dist.yml")

# Random number generator
RNGkind("L'Ecuyer-CMRG")
base_seed <- 15

# Parametrization for the primary and secondary signatures
## Primary
nr_noise <- 0.02
pDEG <- 0.2
p_up=pDEG*dist$alpha
p_down=pDEG*(1-dist$alpha)

p_nb_mol_ctrl <- readRDS("./data/ptt_real.RDS")
p_nb_mol_ctrl <- p_nb_mol_ctrl[!is.na(p_nb_mol_ctrl)]
gene_names <- names(p_nb_mol_ctrl)
nb_genes = length(gene_names)

ind_param <- list(list(p_up= p_up,p_down=p_down,prop_sm_up=dist$prop_sm_up,
                       prop_sm_down=dist$prop_sm_down,
                       qf_vect_up=dist$qf_vect_up,qf_vect_down=dist$qf_vect_down,
                       qf_vect_nr=dist$qf_vect_nr,up_means= dist$up_means,
                       nr_means=dist$nr_means,down_means= dist$down_means))


# Parametrization for the replicates interconnected signatures
cs_inter <- 0.9
nr_noise_inter <- 0.02
submod_transition_inter <- "cop"
copula_submod_inter <- "Frank"
rho_submod_inter <- 0.9
eps_submod_inter <- 1e-3
optim_method_submod_inter <- "Brent"
proba_transition_inter <- "cop"
copula_prob_inter <- "Frank"
theta_prob_inter <- 10
nbins_prob_inter <- 1e3
inter_param = list(mod_transition="default",
                   submod_transition=submod_transition_inter,
                   qf_vect_up=dist$qf_vect_up,
                   qf_vect_down=dist$qf_vect_down,
                   qf_vect_nr=dist$qf_vect_nr,
                   connectivity_score= cs_inter,
                   transition_mat=NA,nr_noise= nr_noise_inter,
                   prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                   proba_transition=proba_transition_inter, 
                   rho_submod= rho_submod_inter, copula_submod = copula_submod_inter,
                   eps_submod = eps_submod_inter,
                   nbins_prob=nbins_prob_inter,
                   optim_method_submod = optim_method_submod_inter,
                   copula_prob = copula_prob_inter,
                   theta_prob = theta_prob_inter)

# Seeds 
set.seed(base_seed)
seeds <- sample(seq(1,5000),7,replace = F)
ncpus <- 7
sim_data <- simulation_pipeline(ncpus=ncpus, nb_genes=nb_genes, 
                                ind_param=ind_param,
                                gene_names=gene_names,
                                dep_param = NULL,
                                p_nb_mol=p_nb_mol_ctrl,
                                inter_param = inter_param,
                                seeds=seeds,return_res = TRUE,
                                save_DEA=NULL,save_lfc_cosimu = NULL,correct_BE = T
)

deseq_simul <- sim_data$deseq_DEA[[1]]
deseq_real <- read.csv("./data/real_lfc_data.csv")

lfc_real <- deseq_real$log2FoldChange
lfc_real <-deseq_real$log2FoldChange[-log10(deseq_real$padj)<lim_padj]
lfc_real <- lfc_real[!is.na(lfc_real)]

lfc_simul <- deseq_simul$log2FoldChange
lfc_simul <-deseq_simul$log2FoldChange[-log10(deseq_simul$padj)<lim_padj]
lfc_simul <- lfc_simul[!is.na(lfc_simul)]

h_real <- hist(lfc_real,breaks = seq(-8,8,0.01),freq = F)
h_sim <- hist(lfc_simul,breaks = seq(-8,8,0.01),freq = F)

kl3 <- KL(rbind(h_sim$counts/sum(h_sim$counts),h_real$counts/sum(h_real$counts)))
# kullback-leibler : 0.1855657
#RATIO 
kl3/mean_kl
# kullback-leibler : 1.519305 

###### KL divergence between two distinct real distributions conditions ########
###### from the same experiment  ######
deseq_real <- read.csv("./data/real_lfc_data.csv")
deseq_real2 <- read.csv("./data/alt_real_lfc_data.csv")

lfc_real <- deseq_real$log2FoldChange
lfc_real <-deseq_real$log2FoldChange[-log10(deseq_real$padj)<lim_padj]
lfc_real <- lfc_real[!is.na(lfc_real)]

lfc_real2 <- deseq_real2$log2FoldChange
lfc_real2 <-deseq_real2$log2FoldChange[-log10(deseq_real2$padj)<lim_padj]
lfc_real2 <- lfc_real2[!is.na(lfc_real2)]

h_real <- hist(lfc_real,breaks = seq(-7.5,6.1,0.01),freq = F)
h_real2 <- hist(lfc_real2,breaks = seq(-7.5,6.1,0.01),freq = F)

kl4 <- KL(rbind(h_real2$counts/sum(h_real2$counts),h_real$counts/sum(h_real$counts)))

#RATIO 
kl4/mean_kl
# kullback-leibler : 5.323613 


# VOLCANO PLOT
volcano_plot <- function(res, x = 'log2FoldChange',y = 'padj', title=NULL){
  df <- as.data.frame(res)
  EnhancedVolcano::EnhancedVolcano(df,
                                   lab = NA,
                                   x = x,
                                   y = y, 
                                   pCutoff = 0.05,
                                   FCcutoff = 1.5,
                                   title = title,
                                   subtitle = NULL,
                                   caption = NULL
  )
}

lim_lfc <- 8

VP_R <- volcano_plot(deseq_real,title = "Real")  + ylim(0,lim_padj) + xlim(-lim_lfc,lim_lfc)

VP_R <- VP_R +
  theme(axis.ticks.length=unit(0.5, "cm"))+ 
  geom_rug(sides = "bl", alpha=0.1, col="gold3", outside = TRUE,
           length = unit(0.5, "cm")) +xlim(-lim_lfc,lim_lfc) 
#VP_R

VP_S <- volcano_plot(deseq_simul, title= "S. with control base") +ylim(0,lim_padj) + xlim(-lim_lfc,lim_lfc) 

VP_S <- VP_S +
  theme(axis.ticks.length=unit(0.5, "cm"))+ 
  geom_rug(sides = "bl", alpha=0.1, col="gold3", outside = TRUE,
           length = unit(0.5, "cm"))+xlim(-lim_lfc,lim_lfc) 

VP_S2 <- volcano_plot(deseq_simul_non_corrected, title = "S. without control base") +ylim(0,lim_padj) + xlim(-lim_lfc,lim_lfc) 

VP_S2 <- VP_S2 +
  theme(axis.ticks.length=unit(0.5, "cm"))+ 
  geom_rug(sides = "bl", alpha=0.1, col="gold3", outside = TRUE,
           length = unit(0.5, "cm"))+xlim(-lim_lfc,lim_lfc) 
#VP_S2
ggarrange(VP_R,VP_S2,VP_S, ncol=3,common.legend = T)



# HISTOGRAM VISUALIZATION
breaks_padj <- seq(0,lim_padj,1)
breaks_lfc <- seq(-lim_lfc,lim_lfc,0.1)
ss_df_simul <- data.frame(log2FoldChange=lfc_simul,data= "s.replicates")
ss_df_rel <- data.frame(log2FoldChange=lfc_real, data= "real")
cosimu_lfc <- data.frame(log2FoldChange = sim_data$list_ind_base[[1]]$get_lfc_vect(), data= "s.single")
df <- rbind(cosimu_lfc,ss_df_simul,ss_df_rel)
df$data <- factor(df$data, levels = c("s.single","s.replicates","real"))
hist_overlap <- ggplot(df,aes(x=log2FoldChange,fill=data))+
  geom_histogram(aes(y=after_stat(density)),position = "identity",breaks = breaks_lfc,alpha=0.6)+
  scale_fill_manual(values = c("#00246b","#74ADD1","#B0B09D"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "top",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = c(0.8,0.8),
        legend.direction = "vertical",
        strip.text = element_text(size = 20),
        legend.title = element_blank())+
  xlim(-4,4)
hist_overlap
