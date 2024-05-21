library(prg)
library(reticulate)
library(purrr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(parallel)
library(dplyr)


use_python("~/miniconda3/bin/python") #need to be update
sklearn<-import(module = "sklearn")
auc<-sklearn$metrics$auc
average_precision_score <- sklearn$metrics$average_precision_score

# FUNCTIONS
AVG_PS <- function(expected_score, estimated_score, N, decreasing=T, interpolated=FALSE){
  topN <- order(expected_score,decreasing = decreasing)[1:N]
  est_rank <- order(estimated_score,decreasing = decreasing)
  if(!interpolated){
    TP <- cumsum(est_rank %in% topN)
    pr <- TP/seq(1,length(est_rank))
    rel <- est_rank %in% topN
    avg_pr <- sum(map_dbl(seq(1,length(pr)),~pr[.x]*rel[.x]/N))
    return(avg_pr)
  }else{
    TP <- cumsum(est_rank %in% topN)
    recall <- TP/N
    thresh <- unique(recall) #thresholds for the recall 
    pr <- TP/seq(1,length(est_rank))
    pr_inter <- map_dbl(seq(0,1,0.001),~max(pr[which(recall>=.x)]))
    return(pr_inter)
  }
}

path_scores <- "./paper_scripts/RDS_asymmetrical/scores/" #"./paper_scripts/RDS_symmetrical/scores/"
df_params <- readRDS("./paper_scripts/RDS_asymmetrical/df_param.RDS") #"./paper_scripts/RDS_symmetrical/df_param.RDS"
CSS <- readRDS(paste0(path_scores,"CSS.rds"))
CMAP2 <- readRDS(paste0(path_scores,"CMAP2.rds"))
CMAP1 <- readRDS(paste0(path_scores,"CMAP1.rds"))
XSum <- readRDS(paste0(path_scores,"XSum_NA.rds"))
XCos<- readRDS(paste0(path_scores,"XCos_NA.rds"))
XCor_spearman <- readRDS(paste0(path_scores,"XCor_spearman_NA.rds"))
XCor_pearson <- readRDS(paste0(path_scores,"XCor_pearson_NA.rds"))


names(CSS) <- names(CMAP2) <- names(CMAP1) <- names(XSum) <- names(XCos) <-
  names(XCor_pearson) <- names(XCor_spearman) <- unique(df_params$ind_id)

ncpus <- 7
color_pal <- c("#6c0000","#F46D43","#00246b","#B50003","#74ADD1","#B0B09D","#0057c4")
names(color_pal) <- c("cmap1",
                      "cmap2",
                      "css",
                      "xcos",
                      "xpearson",
                      "xspearman",
                      "xsum")
shape = NULL
scores <- NULL

generate_res <- function(df_params,noise_val,method, decreasing, N){
  list_res <- mclapply(unique(df_params$rep_id),function(i){
    ids <- unique(df_params$ind_id[df_params$rep_id==i])
    cmap1 <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = CMAP1[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    cmap2 <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = CMAP2[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    css <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = CSS[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    xsum <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = XSum[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    xcos <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = XCos[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    xpearson <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = XCor_pearson[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    xspearman <- bind_cols(lapply(ids, function(s){
      map_dbl(N,
              function(n){
                tmp <- df_params[df_params$ind_id==s,]
                sel <- which(tmp$nr_noise %in% noise_val)
                method(expected_score = tmp$cs[sel],
                       estimated_score = XCor_spearman[[s]][sel],N = n,decreasing = decreasing)
              })
    }))
    data <- rbind(cmap1,
                  cmap2,
                  css,
                  xsum,
                  xcos,
                  xpearson,
                  xspearman)
    nb_col <- ncol(data)
    data$mean <- rowMeans(data)
    data$sd <- map_dbl(seq(1,nrow(data)), ~sd(data[.x,1:nb_col]))
    data$max_ci <- data$mean + 1.96*data$sd/sqrt(nb_col)
    data$min_ci <- data$mean - 1.96*data$sd/sqrt(nb_col) 
    #data$min_ci[data$min_ci<0] <- 0
    data$N <- rep(N,nrow(data)/length(N))
    data$estimator <- rep(c("cmap1",
                            "cmap2",
                            "css",
                            "xsum",
                            "xcos",
                            "xpearson",
                            "xspearman"),each=length(N))
    p <- ggplot(data, aes(x=as.factor(N),y=mean,fill=estimator))+
      scale_fill_manual(values=color_pal) +
      geom_bar(stat="identity", color="black", position=position_dodge())+
      #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), stat = "identity",position=position_dodge(1))+
      geom_errorbar(aes(ymin=min_ci, ymax=max_ci), stat = "identity",position=position_dodge())+
      labs(color="estimator")+
      xlab("Top N") + ylab("Avg. Precision Score")+
      theme_classic() + 
      theme(panel.background = element_rect(fill = "#FFFFFF"),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.position = "top")+ylim(0,1) +
      labs(fill="")+scale_fill_manual(values=color_pal)
    
  })
  return(list_res)
}

df_params$rep_id <- 1
df_params$set <- 1
N <- c(1,3,5,10,20,50)
decreasing <- TRUE

# AVG_PS
method <- AVG_PS
list_res <- generate_res(df_params = df_params,noise_val = c(0.02), method= method, decreasing = decreasing,
                         N = N)
list_res[[1]]

