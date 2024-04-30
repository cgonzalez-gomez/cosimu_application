#### ATTENTION : the paths are noy updated in this script ####

library(purrr)

# path <- "./RDS_asymmetrical/"
path <- "./RDS_symmetrical/"

df_param <- readRDS(paste0(path,"df_param.RDS"))
df_padj_queries <- map_dfc(unique(df_param$ind_base), function(i){
  deseq_data <- readRDS(paste0(path,"Q_",i,"/deseq2.RDS"))
  df <- data.frame(deseq_data[[1]]$padj)
  colnames(df) <- paste0("ind_",i)
  return(df)
})

df_lfc_queries <- map_dfc(unique(df_param$ind_base), function(i){
  deseq_data <- readRDS(paste0(path,"Q_",i,"/deseq2.RDS"))
  df <- data.frame(deseq_data[[1]]$log2FoldChange)
  colnames(df) <- paste0("ind_",i)
  return(df)
})

df_lfc_profiles <- map_dfc(unique(df_param$ind_base), function(i){
  deseq_data <- readRDS(paste0(path,"Q_",i,"/deseq2.RDS"))
  df <- map_dfc(seq(2,length(deseq_data)), ~deseq_data[[.x]]$log2FoldChange)
  return(df)
})
colnames(df_lfc_profiles) <- paste0("dep_",seq(1,ncol(df_lfc_profiles)))

saveRDS(df_lfc_profiles,paste0(path,"mat_profiles.RDS"))
saveRDS(df_lfc_queries,paste0(path,"mat_queries.RDS"))
saveRDS(df_padj_queries,paste0(path,"mat_queries_padj.RDS"))
