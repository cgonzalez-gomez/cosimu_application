# IMPORTANT NOTE: the script is not updated with the last version of the simulation pipeline #

library(readr)
library(stringr)
library(DESeq2)
library(purrr)
library(GEOquery)

path_data <- "./downloaded_real_data/" # need UPDATE
# Data Download
# GSE181472
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE181472&format=file&file=GSE181472%5F2020%5F1005%5FfeatureCounts%5Fgenes%2Etsv%2Egz",
              destfile = paste0(path_data,"GSE181472_2020_1005_featureCounts_genes.tsv"), mode = "wb")
download.file(url ="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181472/soft/GSE181472_family.soft.gz",
              destfile = paste0(path_data,"GSE181472_family.soft"), mode = "wb")
# GSE182024
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE182024&format=file&file=GSE182024%5FRaw%5Fgene%5Fcounts%5Fmatrix%2Etxt%2Egz",
              destfile = paste0(path_data,"GSE182024_Raw_gene_counts_matrix.txt"), mode = "wb")
download.file(url ="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182024/soft/GSE182024_family.soft.gz",
              destfile = paste0(path_data,"GSE181472_family.soft"), mode = "wb")
# GSE185453
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE185453&format=file&file=GSE185453%5FRNA%5Fraw%5Fcounts%5Fmatrix%2Etxt%2Egz",
              destfile = paste0(path_data,"GSE185453_RNA_raw_counts_matrix.txt"), mode = "wb")
download.file(url ="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185453/soft/GSE185453_family.soft.gz",
              destfile = paste0(path_data,"GSE185453_family.soft"), mode = "wb")


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
data_GSE181472 <- getGEO(filename = paste0(path_data,"GSE181472_family.soft"))

sample_id <- data_GSE181472@header$sample_id
description <- map_chr(sample_id,~data_GSE181472@gsms[[.x]]@header$description)

cts_df <- read.table(file = paste0(path_data,'GSE181472_2020_1005_featureCounts_genes.tsv'), sep = '\t', header = TRUE)

colnames(cts_df)<-str_remove(colnames(cts_df),"X")
cts_df <- cts_df[,description]

cts_df<-cts_df[rowSums(cts_df)!=0,]
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
# 2 replicates only

data_GSE185453 <- getGEO(filename = paste0(path_data,"GSE185453_family.soft"))
cts_df <- read.table(file = paste0(path_data,'GSE185453_RNA_raw_counts_matrix.txt'), sep = '\t', header = TRUE)
cts_df<-cts_df[rowSums(cts_df)!=0,]
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

saveRDS(res_GSE181472,paste0(path_data,"res_GSE181472.RDS"))
saveRDS(res_GSE182024,paste0(path_data,"res_GSE182024.RDS"))
saveRDS(res_GSE185453,paste0(path_data,"res_GSE185453.RDS"))
