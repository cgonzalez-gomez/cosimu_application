####### COMPARISON VS EXPERIMENTAL REPLICATES #######
# Correlation analysis, with simulated data (c=1) #

library(mgcv)
library(ggplot2)
library(ggnewscale)

source("./paper_scripts/simulation_pipeline.R")
path_data <- "./downloaded_real_data/" #NEED UPDATE
res_GSE181472 <- readRDS(paste0(path_data,"res_GSE181472.RDS"))
res_GSE182024 <- readRDS(paste0(path_data,"res_GSE182024.RDS"))
res_GSE185453 <- readRDS(paste0(path_data,"res_GSE185453.RDS")) #two replicates

submod_transition_inter <- "deterministic"#"cop"
submod_transition <- "deterministic" #"cop"
proba_transition_inter <- "deterministic"
proba_transition <-"deterministic"
nr_noise_inter <- 0.02
nr_noise <- 0.02
size_factor <- 0.05

### SIMULATED DATA 
dist <- cosimu::load_dist("./data/real_symmetrical_dist.yml")
RNGkind("L'Ecuyer-CMRG")

alpha <- 0.5

cs_inter <- 1

cs <- 1


p_nb_mol_ctrl <- readRDS("./data/ptt_real.RDS")
p_nb_mol_ctrl <- p_nb_mol_ctrl[!is.na(p_nb_mol_ctrl)]
gene_names <- names(p_nb_mol_ctrl)
nb_genes = length(gene_names)

nb_sign <- nrow(res_GSE181472)+nrow(res_GSE182024)+nrow(res_GSE185453)
# 
# x <- runif(min=0,max=0.4,n = nb_sign)
pDEG_vect <- rep(c(0.01,0.05,0.15,0.25,0.35), each=4)
inter_param = list(mod_transition="default",
                   submod_transition=submod_transition_inter,
                   qf_vect_up=dist$qf_vect_up,
                   qf_vect_down=dist$qf_vect_down,
                   qf_vect_nr=dist$qf_vect_nr,
                   connectivity_score= cs_inter,
                   transition_mat=NA,nr_noise= nr_noise_inter,
                   prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                   proba_transition=proba_transition_inter)

dep_param = list(list(mod_transition="default",
                      submod_transition=submod_transition,
                      qf_vect_up=dist$qf_vect_up,
                      qf_vect_down=dist$qf_vect_down,
                      qf_vect_nr=dist$qf_vect_nr,
                      connectivity_score= cs,
                      transition_mat=NULL,nr_noise= nr_noise,
                      prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                      proba_transition=proba_transition,up_means= dist$up_means,
                      nr_means=dist$nr_means,down_means= dist$down_means))
ncpus <- 7

set.seed(100)
tictoc::tic()
res_simulated <- map_dfr(pDEG_vect,function(pDEG){
  seeds <-sample(seq(1,5000),7,replace = F)
  ind_param = list(list(p_up= pDEG*alpha,p_down=pDEG*(1-alpha),prop_sm_up=dist$prop_sm_up,prop_sm_down=dist$prop_sm_down,
                        qf_vect_up=dist$qf_vect_up,qf_vect_down=dist$qf_vect_down,qf_vect_nr=dist$qf_vect_nr,up_means= dist$up_means,
                        nr_means=dist$nr_means,down_means= dist$down_means))
  
  res <- simulation_pipeline(ncpus=ncpus, nb_genes=nb_genes, 
                             ind_param=ind_param,
                             gene_names=gene_names,
                             dep_param = dep_param,
                             p_nb_mol=p_nb_mol_ctrl,
                             inter_param = inter_param,
                             seeds=seeds,
                             save_lfc_cosimu = NULL,
                             save_DEA = NULL,size_factor = size_factor,
                             ncpus_secondary = 1, paired= F,
                             return_res = TRUE)
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
tictoc::toc()


saveRDS(res_simulated, paste0(path_data,"res_simulated_grp.RDS"))
class <- rep(c("real","simulated"),c(nb_sign,
                                     nrow(res_simulated)))
set <- c(rep("GSE181472",nrow(res_GSE181472)),
         rep("GSE182024",nrow(res_GSE182024)),
         rep("GSE185453",nrow(res_GSE185453)), 
         rep("simulation", nrow(res_simulated) ))
glob_res <- rbind(res_GSE181472,
                  res_GSE182024,
                  res_GSE185453,
                  res_simulated)
glob_res<-cbind(glob_res,set=as.factor(set),class=as.factor(class))
glob_res$meanDEG.lfc.padj <- rowMeans(glob_res[,c("p.DEG1","p.DEG2")])

#PLOT
color_pal <- c("#F46D43","#B50003","#6c0000","#74ADD1")
line_pal<- c("#B50003","#74ADD1")

g <- ggplot(glob_res,aes(x=meanDEG.lfc.padj,y=cor))+
  geom_point(aes(color=set),size=3)+
  scale_shape_manual(values = c(rep(16,3),rep(17,5)))+
  scale_color_manual(values=color_pal, name = "Data Source")+
  new_scale_colour() +
  geom_smooth(aes(group=class, color=class),method = "gam",
              size=0.5, formula = y ~ s(x, bs = "cs"),show.legend = F,se=F)+
  scale_color_manual(values=line_pal)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "right",
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 25))+
  xlab("Mean proportion of DEG")+
  ylab("Pearson's correlation")+
  labs(shape=NULL, color=NULL)
g

# Kruskall
glob_res$DEG_class <- map_chr(glob_res$meanDEG.lfc.padj,function(x){
  if(x < 0.02){
    return("low")
  }else if(x>0.04){
    return("high")
  }else{
    return("medium")
  }
})

kruskal.test(cor ~ DEG_class, data=glob_res[glob_res$class=="real",])
kruskal.test(cor ~ DEG_class, data=glob_res[glob_res$class=="simulated",])
