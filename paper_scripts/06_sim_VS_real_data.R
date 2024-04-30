####### COMPARISON VS EXPERIMENTAL REPLICATES #######
# Correlation analysis, with simulated data (c=1) #
# IMPORTANT NOTE: the script is not updated with the last version of the simulation pipeline #

library(readr)
library(stringr)
library(DESeq2)
library(purrr)
library(GEOquery)
library(ggplot2)


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


compare <- function(df,untreat_cols, treat_cols){
  stopifnot(length(untreat_cols)==length(treat_cols))
  N <- length(untreat_cols)
  cts1 <- df[,c(untreat_cols[1:(N/2)],treat_cols[1:(N/2)])]
  cts2 <- df[,c(untreat_cols[(N/2+1):N],treat_cols[(N/2+1):N])]
  num_reps <- N/2
  design = "~ condition"
  
  coldata1 <- data.frame(condition = rep(c("untreated", "treated"), each = num_reps))
  rownames(coldata1) <- colnames(cts1)
  
  coldata2 <- data.frame(condition = rep(c("untreated", "treated"), each = num_reps))
  rownames(coldata2) <- colnames(cts2)
  
  res1 <- deseq_analysis(cts1,coldata1,design)
  res2 <- deseq_analysis(cts2,coldata2,design)
  
  return(data.frame(cor=cor(res1$log2FoldChange,res2$log2FoldChange,use = "na.or.complete"),
                    p.padj1=sum(res1$padj<0.05,na.rm=T)/nrow(res1),
                    p.padj2=sum(res2$padj<0.05,na.rm=T)/nrow(res2),
                    p.lfc1=sum(abs(res1$log2FoldChange) > 1.5,na.rm=T)/nrow(res1),
                    p.lfc2=sum(abs(res2$log2FoldChange) > 1.5,na.rm=T)/nrow(res2),
                    p.DEG1 = sum(abs(res1$log2FoldChange) > 1.5 & res1$padj<0.05,na.rm=T)/nrow(res1),
                    p.DEG2 = sum(abs(res2$log2FoldChange) > 1.5 & res2$padj<0.05,na.rm=T)/nrow(res2)
                    
  ))
}

###### REAL DATA ######
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181472
# Effects of nintedanib on unstimulated and stimulated macrophages 
# from human monocytes in single culture and in co-culture with fibroblasts.

data_GSE181472 <- getGEO(filename = "/exports/home/cgonzalez/cosimu_paper/new_version_12_2022/real_data/GSE181472_family.soft")
#View(data_GSE181472)

sample_id <- data_GSE181472@header$sample_id
description <- map_chr(sample_id,~data_GSE181472@gsms[[.x]]@header$description)

cts_df <- read.table(file = '/exports/home/cgonzalez/cosimu_paper/new_version_12_2022/real_data/GSE181472_2020_1005_featureCounts_genes.tsv', sep = '\t', header = TRUE)
colnames(cts_df)<-str_remove(colnames(cts_df),"X")
cts_df <- cts_df[,description]

sel <- order(rowSums(cts_df),decreasing = T,na.last = T)
cts_df<-cts_df[sel[1:10000],]
characteristics <- map_dfr(sample_id,function(x) {
  row <- str_remove(data_GSE181472@gsms[[x]]@header$characteristics_ch1,".*: ")
  names(row) <- c("treatment", "cult_cond")
  NTB <- "NTB" %in% str_split(row["treatment"],"_")[[1]]
  TNFa <- "TNFa" %in% str_split(row["treatment"],"_")[[1]]
  return(c(row,NTB=NTB,TNFa=TNFa))
})

characteristics$description <- description
characteristics$cult_cond <- as.factor(characteristics$cult_cond)

combinations <- unique(characteristics[,c("cult_cond","TNFa")]) 

res_GSE181472 <- map_dfr(1:nrow(combinations),function(i){# Comparison of NTB (at comparable level of TNF)
  cond <- combinations[i,]
  ctrl_cols <- characteristics$description[which(characteristics$NTB == FALSE &
                                                   characteristics$cult_cond == cond$cult_cond &
                                                   characteristics$TNFa ==  cond$TNFa)]
  treat_cols <- characteristics$description[which(characteristics$NTB == TRUE &
                                                    characteristics$cult_cond == cond$cult_cond &
                                                    characteristics$TNFa ==  cond$TNFa)]
  compare(cts_df,ctrl_cols,treat_cols)
})

res_GSE181472 <- rbind(res_GSE181472, 
                       map_dfr(levels(characteristics$cult_cond),function(l){ # Comparison TNF AND NTB
                         ctrl_cols <- characteristics$description[which(characteristics$NTB == FALSE &
                                                                          characteristics$cult_cond == l &
                                                                          characteristics$TNFa ==  FALSE)]
                         treat_cols <- characteristics$description[which(characteristics$NTB == TRUE &
                                                                           characteristics$cult_cond == l &
                                                                           characteristics$TNFa ==  TRUE)]
                         compare(cts_df,ctrl_cols,treat_cols)
                       }),
                       map_dfr(levels(characteristics$cult_cond),function(l){ # Comparison of TNF (without NTB)
                         ctrl_cols <- characteristics$description[which(characteristics$NTB == FALSE &
                                                                          characteristics$cult_cond == l &
                                                                          characteristics$TNFa ==  FALSE)]
                         treat_cols <- characteristics$description[which(characteristics$NTB == FALSE &
                                                                           characteristics$cult_cond == l &
                                                                           characteristics$TNFa ==  TRUE)]
                         compare(cts_df,ctrl_cols,treat_cols)
                       })
)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182024
data_GSE182024 <- getGEO(filename = "/exports/home/cgonzalez/cosimu_paper/new_version_12_2022/real_data/GSE182024_family.soft")

sample_id <- data_GSE182024@header$sample_id
description <- map_chr(sample_id,~data_GSE182024@gsms[[.x]]@header$title)

cts_df <- read.table(file = '/exports/home/cgonzalez/cosimu_paper/new_version_12_2022/real_data/GSE182024_Raw_gene_counts_matrix.txt', sep = '\t', header = TRUE)
cts_df <- cts_df[,description]


characteristics <- map_dfr(description,function(x) {
  row <- str_split(x,"_")[[1]]
  ctrl <- "Control" %in% row
  TNFa <- "TNFa" %in% row
  return(c(ctrl=ctrl,IL4_IL13 = (!ctrl & !TNFa),TNFa=TNFa))
})


characteristics$description <- description

res_GSE182024 <- map_dfr(c("IL4_IL13","TNFa"),function(x){
  ctrl_cols <- characteristics$description[which(characteristics$ctrl==TRUE)]
  treat_cols <- characteristics$description[which(characteristics[,x] == TRUE)]
  compare(cts_df,ctrl_cols,treat_cols)
})

res_GSE182024<- rbind(res_GSE182024,compare(cts_df,
                                            untreat_cols = characteristics$description[which(characteristics$IL4_IL13 == TRUE)],
                                            treat_cols = characteristics$description[which(characteristics$TNFa == TRUE)]
))


# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185453 
# uniquement 2 réplicats 

data_GSE185453 <- getGEO(filename = "/exports/home/cgonzalez/cosimu_paper/new_version_12_2022/real_data/GSE185453_family.soft")
cts_df <- read.table(file = '/exports/home/cgonzalez/cosimu_paper/new_version_12_2022/real_data/GSE185453_RNA_raw_counts_matrix.txt', sep = '\t', header = TRUE)

characteristics <- map_dfr(colnames(cts_df), function(x){
  row <- str_split(x,"[_.]")[[1]]
  names(row) <- c("type","region","treat1","treat2","rep")
  return(row)
}) 
characteristics$description <- colnames(cts_df)
res_GSE185453<- rbind(compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "CFC")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "Context")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "Veh" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "Veh" & characteristics$treat2 == "Context")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "Veh" & characteristics$treat2 == "CFC")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "Veh" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "Veh" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "Context")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "Veh" & characteristics$treat2 == "CFC")],
                              treat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "Veh" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "hippocampus" & characteristics$treat1 == "Veh" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "Veh" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "Context")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "Veh" & characteristics$treat2 == "CFC")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "HDACi" & characteristics$treat2 == "CFC")]),
                      compare(cts_df,
                              untreat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "Veh" & characteristics$treat2 == "Context")],
                              treat_cols = characteristics$description[which(characteristics$region == "striatum" & characteristics$treat1 == "Veh" & characteristics$treat2 == "CFC")])
)

saveRDS(res_GSE181472,"~/fig_corr1/res_GSE181472.RDS")
saveRDS(res_GSE182024,"~/fig_corr1/res_GSE182024.RDS")
saveRDS(res_GSE185453,"~/fig_corr1/res_GSE185453.RDS")

#########################
res_GSE181472 <- readRDS("~/fig_corr1/res_GSE181472.RDS")
res_GSE182024 <- readRDS("~/fig_corr1/res_GSE182024.RDS")
res_GSE185453 <- readRDS("~/fig_corr1/res_GSE185453.RDS") #two replicates

load("~/fig2/mean_modes_dist.Rdata")
RNGkind("L'Ecuyer-CMRG")
### SIMULATED DATA 

mean_gamma <- function(shape,scale) shape*scale
source("/exports/analysis/PROJECTS/SIGNATURA_CONNECT/formal_var_models/formal_sim_pipeline.R")

alpha <- 0.5

cs_inter <- 1
nr_noise_inter <- 0.02
submod_transition_inter <- "cop" #"deterministic"#
copula_submod_inter <- "Frank"
rho_submod_inter <- 0.9
eps_submod_inter <- 1e-3
optim_method_submod_inter <- "Brent"
perc_transition_inter <- "stochastic" #"deterministic"#
beta_error_inter <- 5
inter_param = list(p_up=NA,p_down=NA, mod_transition="nocop",
                   submod_transition=submod_transition_inter,
                   qf_vect_up=dist$qf_vect_up,
                   qf_vect_down=dist$qf_vect_down,qf_vect_nr=dist$qf_vect_nr,expected_cs= cs_inter,
                   transition_mat=NULL,nr_noise= nr_noise_inter,
                   prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                   perc_transition=perc_transition_inter , beta_error =beta_error_inter,
                   rho_submod= rho_submod_inter, copula_submod = copula_submod_inter,
                   eps_submod = eps_submod_inter,
                   optim_method_submod = optim_method_submod_inter)


sigma_bio <- 0.5
sigma_poly <- 3
num_reps <- 3
correct_BE <- TRUE
num_reads=4300000



cs <- 1
nr_noise <- 0.02
submod_transition <- "cop"
copula_submod <- "Frank"
rho_submod <- 0.9
eps_submod <- 1e-3
optim_method_submod <- "Brent"
perc_transition <- "stochastic"
beta_error <- 5
dep_param = list(list(p_up=NA,p_down=NA, mod_transition="nocop",
                      submod_transition=submod_transition,
                      qf_vect_up=dist$qf_vect_up,
                      qf_vect_down=dist$qf_vect_down,qf_vect_nr=dist$qf_vect_nr,expected_cs= cs,
                      transition_mat=NULL,nr_noise= nr_noise,
                      prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                      perc_transition=perc_transition, beta_error =beta_error,
                      rho_submod= rho_submod, copula_submod = copula_submod,
                      eps_submod = eps_submod,
                      optim_method_submod = optim_method_submod,up_means= dist$up_means,
                      nr_means=dist$nr_means,down_means= dist$down_means))

ncpus <- 1
# gene_names = readRDS("/exports/analysis/PROJECTS/SIGNATURA_CONNECT/annotations/gene_names.RDS")
# p_nb_mol_base <- readRDS("/exports/analysis/PROJECTS/SIGNATURA_CONNECT/annotations/ptt_baseline_genes_D10_04_2023.RDS")[gene_names]
# n <-10000
# p_nb_mol_base <- sort(p_nb_mol_base,decreasing = T)[1:n]
# gene_names <- names(p_nb_mol_base)
# nb_genes = length(gene_names)


p_nb_mol_base <- readRDS("~/fig2/ptt_real.RDS")
p_nb_mol_base <- p_nb_mol_base[!is.na(p_nb_mol_base)]
n <- 7000
p_nb_mol_base <- sort(p_nb_mol_base,decreasing = T)[1:n]
gene_names <- names(p_nb_mol_base)
nb_genes = length(gene_names)


NUMBER_REPETITIONS <- 4
pDEG_vect <- c(0,0.01,0.05,0.1,0.15)
paired = TRUE
set.seed(100)
tictoc::tic()
res_simulated <- map(pDEG_vect,function(pDEG){
  map_dfr(seq(1,NUMBER_REPETITIONS),function(i){
    seeds <-sample(seq(1,1000),6)
    ind_param = list(list(p_up= pDEG*alpha,p_down=pDEG*(1-alpha),prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                          qf_vect_up=dist$qf_vect_up,qf_vect_down=dist$qf_vect_down,qf_vect_nr=dist$qf_vect_nr,up_means= dist$up_means,
                          nr_means=dist$nr_means,down_means= dist$down_means))
    res <- full_pipe(ncpus=ncpus, nb_genes=nb_genes, ind_param=ind_param,gene_names=gene_names,dep_param = dep_param,
                     variability = "biological",interaction = TRUE, p_nb_mol=p_nb_mol_base, num_reps= num_reps,inter_param = inter_param,
                     df_base_fc = NULL, DE_analysis = "deseq",counts_only = TRUE, ctrl_fc = "ID",
                     return_res = TRUE,get_nb_DEG = F,sigleton_analysis = FALSE,
                     size_factor = sigma_poly,save_lfc_cosimu = NULL,save_counts_pol = NULL,save_DEA = NULL,
                     save_DEG = NULL,save_singleton = NULL,sigma_bio = sigma_bio,correct_BE = correct_BE, seeds=seeds, paired = paired,
                     add_low = F,basemeans_thr = 1)
    res1 <- res$deseq_DEA[[1]]
    res2<-res$deseq_DEA[[2]]
    return(data.frame(cor=cor(res1$log2FoldChange,res2$log2FoldChange,use = "na.or.complete"),
                      p.padj1=sum(res1$padj<0.05,na.rm=T)/nrow(res1),
                      p.padj2=sum(res2$padj<0.05,na.rm=T)/nrow(res2),
                      p.lfc1=sum(abs(res1$log2FoldChange) > 1.5,na.rm=T)/nrow(res1),
                      p.lfc2=sum(abs(res2$log2FoldChange) > 1.5,na.rm=T)/nrow(res2),
                      p.DEG1 = sum(abs(res1$log2FoldChange) > 1.5 & res1$padj<0.05,na.rm=T)/nrow(res1),
                      p.DEG2 = sum(abs(res2$log2FoldChange) > 1.5 & res2$padj<0.05,na.rm=T)/nrow(res2)))
    
  })
})
tictoc::toc()
# 1 pDEG_vect <- c(0,0.01,0.05,0.1) ;NUMBER_REPETITIONS : 8; proba_trans: stochastic ; sm_trans: cop
# 2 pDEG_vect <- c(0,0.01,0.05,0.1,0.2,0.25) ; NUMBER_REPETITIONS: 4; proba_trans: stochastic ; sm_trans: cop
# 3 pDEG_vect <- c(0,0.01,0.05,0.1,0.2,0.25) ; NUMBER_REPETITIONS: 4; proba_trans: deter ; sm_trans: cop
# 4 pDEG_vect <- c(0,0.01,0.05,0.1,0.2,0.25) ; NUMBER_REPETITIONS: 4; proba_trans: deter; sm_trans: deter
# 5 pDEG_vect <- c(0,0.01,0.05,0.1,0.15);  NUMBER_REPETITIONS: 4; proba_trans: deter; sm_trans: deter
#saveRDS(res_simulated,"~/fig_corr1/res_simulated5.RDS")

saveRDS(res_simulated,"~/fig_corr1/res_simulated_paired_7000.RDS")

#### PLOTS

class <- rep(c("real","simulated"),c(nrow(res_GSE181472)
                                     +nrow(res_GSE182024)
                                     +nrow(res_GSE185453)
                                     ,
                                     length(res_simulated)*NUMBER_REPETITIONS))
set <- c(rep("GSE181472",nrow(res_GSE181472)),
         rep("GSE182024",nrow(res_GSE182024)),
         rep("GSE185453",nrow(res_GSE185453)), # 2 reps
         rep(paste0("simul",seq(1,length(pDEG_vect))),each=NUMBER_REPETITIONS)) 
glob_res <- rbind(res_GSE181472,
                  res_GSE182024,
                  res_GSE185453,
                  map_dfr(res_simulated,~.x))

glob_res<-cbind(glob_res,set=as.factor(set),class=as.factor(class))
glob_res$meanDEG.lfc.padj <- rowMeans(glob_res[,c("p.DEG1","p.DEG2")])
# glob_res$meanDEG.padj <- rowMeans(glob_res[,c("p.padj1","p.padj2")])
# glob_res$meanDEG.lfc <- rowMeans(glob_res[,c("p.lfc1","p.lfc2")])
# 
# ggplot(glob_res,aes(x=meanDEG.padj,y=cor,shape=class,fill=set,color=set))+
#   geom_point()+
#   geom_smooth(aes(group=class),method="loess",color="gray",
#               size=0.5,alpha=0.2)
# 
# ggplot(glob_res,aes(x=meanDEG.lfc,y=cor,shape=class,fill=set,color=set))+
#   geom_point()+
#   geom_smooth(aes(group=class),method="gam",color="gray",
#               size=0.5,alpha=0.2)
# 
# ggplot(glob_res,aes(x=meanDEG.lfc.padj,y=cor,shape=class,fill=set,color=set))+
#   geom_point()+
#   geom_smooth(aes(group=class),method="loess",color="gray",
#               size=0.5,alpha=0.2)
color_pal <- c("#F46D43","#B50003","#6c0000","#0D0556","#00246b","#0057c4","#74ADD1","#6DE2E6")
ggplot(glob_res,aes(x=meanDEG.lfc.padj,y=cor))+
  geom_point(aes(shape=set,color=set),size=3)+
  scale_shape_manual(values = c(rep(16,3),rep(17,5)))+
  scale_color_manual(values=color_pal)+
  geom_smooth(aes(group=class),method = "gam",color="gray",
              size=0.5,alpha=0.2, formula = y ~ s(x, bs = "cs"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "right",
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 25))+
  xlab("Mean proportion of DEG")+
  ylab("Pearson's correlation")+
  labs(shape=NULL, color=NULL)

