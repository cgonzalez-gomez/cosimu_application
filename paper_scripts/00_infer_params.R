########### INFER PARAMETRIZATION FROM REAL DATA##########

real_lfc_data <- read.csv("./data/real_lfc_data.csv")
cosimu::infer_dist(lfc_vector = real_lfc_data$log2FoldChange,save_yaml = "./data/real_asymmetrical_dist.yml")
cosimu::infer_dist(lfc_vector = real_lfc_data$log2FoldChange,save_yaml = "./data/real_symmetrical_dist.yml",sym_mode = TRUE)


########### PROPORTIONAL TO TRUE NUMBER OF TRANSCRIPTS #########
library(AnnotationHub)
cts_mat<-read.table(file = '/exports/analysis/PROJECTS/SIGNATURA_CONNECT/formal_var_models/var_degree/real_data/GSE185985_rna-seq_count_matrix.txt',
                    sep = '\t', row.names = 1,header = T)
cts_mat <- cts_mat[,1:9]
annotation <- map_dfr(seq(1,nrow(cts_mat)),function(i){
  res <- str_split(rownames(cts_mat)[i],":")[[1]]
  names(res) <- c("gene_id_version","tx_id_version","Region")
  gene_id <- str_split(res["gene_id_version"],"\\.")[[1]][1]
  tx_id <- str_split(res["tx_id_version"],"\\.")[[1]][1]
  return(c(res,gene_id=gene_id,tx_id=tx_id))
})

hub <- AnnotationHub()
query(hub, c("homo sapiens","ensdb"))

ensdb <- hub[["AH109606"]]
require("ensembldb")

##Genes
gns <- genes(ensdb)
gns <- gns[gns$gene_id %in% annotation$gene_id,]
gnwid <- setNames(width(gns), names(gns))
gnwid <- gnwid[annotation$gene_id]


#fpk norm 
# fpk_fact <- 1e+03/gnwid
# x <- rowMeans(cts_mat[,1:3]*fpk_fact, na.rm=T)

# TPM intra norm 
tpm_fact<- sweep(matrix(1e+06/gnwid,ncol=3,nrow=length(gnwid)),2,colSums(cts_mat[,seq(1,3)]/gnwid,na.rm = T),FUN = "/")

norm_mat <- cts_mat[,1:3]*tpm_fact

ptt_basemean <- rowMeans(norm_mat,na.rm = T)
names(ptt_basemean) <- annotation$gene_id
saveRDS(ptt_basemean,file="~/fig2/ptt_real.RDS")
saveRDS(gnwid, file = "~/fig2/gene_width.RDS")
