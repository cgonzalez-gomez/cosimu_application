library(reshape2)
library(ggpubr)
library(purrr)
# parametrization 
df_param <- readRDS("./RDS_asymmetrical/df_param.RDS")

# MODELS
path_scores <- "./RDS_asymmetrical/scores/" 

CSS <- readRDS(paste0(path_scores,"CSS.rds"))
CMAP2 <- readRDS(paste0(path_scores,"CMAP2.rds"))
CMAP1 <- readRDS(paste0(path_scores,"CMAP1.rds"))
XSum_NA <- readRDS(paste0(path_scores,"XSum_NA.rds"))
XCos_NA <- readRDS(paste0(path_scores,"XCos_NA.rds"))
XSpe_NA <- readRDS(paste0(path_scores,"XCor_spearman_NA.rds"))
XPear_NA <- readRDS(paste0(path_scores,"XCor_pearson_NA.rds"))
n <- 20

color_pal <- c("#6c0000","#F46D43","#00246b","#B50003","#74ADD1","#B0B09D","#0057c4")
names(color_pal) <- c("cmap1",
                      "cmap2",
                      "css",
                      "xcos",
                      "xpearson",
                      "xspearman",
                      "xsum")

CS_df <- map_dfr(seq(1,n), function(i){
  data.frame(cmap1 = CMAP1[[i]], cmap2 = CMAP2[[i]],
             css = CSS[[i]],
             xsum = XSum_NA[[i]],
             xcos = XCos_NA[[i]],
             xpearson = XPear_NA[[i]],
             xspearman = XSpe_NA[[i]]
  )
  
})
CS_df <- melt(CS_df)  
nb_scores <- 7
CS_df$cs_expected<- rep(df_param$cs,times = nb_scores)
colnames(CS_df) <- c("score", "estimate","cs_expected")

scores <- c(css="CSS", cmap2="CMAP2",cmap1="CMAP1", xsum="XSum_NA",xcos="XCos_NA",
            xspearman="XSpe_NA", xpearson="XPear_NA")
scores <- scores[c("cmap1",
                   "cmap2",
                   "css",
                   "xsum",
                   "xcos",
                   "xpearson",
                   "xspearman")]
mean_df <- map_dfr(names(scores),
                   function(s){
                     data <- as.data.frame(do.call(cbind, get(scores[s])))
                     nb_col <- ncol(data)
                     res <- data.frame(mean= rowMeans(data),sd=map_dbl(seq(1,nrow(data)), ~sd(data[.x,1:nb_col])))
                     res$max_ci <- res$mean + 1.96*res$sd/sqrt(nb_col)
                     res$min_ci <- res$mean - 1.96*res$sd/sqrt(nb_col)
                     res$score <- s
                     res$id <- seq(1,nrow(res))
                     return(res)
                   })
mean_df$cs_expected <- unique(df_param$cs)
min_max <- function(vect) (vect - min(vect))/(max(vect)-min(vect))

mean_df_norm <- map_dfr(names(scores),
                        function(s){
                          data <- as.data.frame(do.call(cbind, map(get(scores[s]),~min_max(.x))))
                          nb_col <- ncol(data)
                          res <- data.frame(mean= rowMeans(data),sd=map_dbl(seq(1,nrow(data)), ~sd(data[.x,1:nb_col])))
                          res$max_ci <- res$mean + 1.96*res$sd/sqrt(nb_col)
                          res$min_ci <- res$mean - 1.96*res$sd/sqrt(nb_col)
                          res$score <- s
                          res$id <- seq(1,nrow(res))
                          return(res)
                        })
mean_df_norm$cs_expected <- unique(df_param$cs)
mean_df_norm$score <- factor(mean_df_norm$score, levels = names(scores))


# Bayesian credible interval
final_plt <- ggplot(data=mean_df_norm, mapping=aes(x = cs_expected,y = mean,color = score, fill=score))+
  geom_point(alpha = 1, size = 0.7)+
  geom_smooth(method = "gam",se = TRUE,size =0.5)+
  xlab("Expected score") + ylab("Normalized estimator")+
  theme_classic() +
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        axis.text.x = element_text(size = 11,hjust = 0.6),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = c(0.9,0.25))+
  ylim(0,1)+
  facet_wrap(~score,scales="free_y",nrow = 2,ncol = 4)+
  scale_fill_manual(name="score",values=color_pal)+
  scale_color_manual(name="score",values=color_pal)

final_plt
