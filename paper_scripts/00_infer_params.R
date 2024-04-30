########### INFER PARAMETRIZATION FROM REAL DATA##########

real_lfc_data <- read.csv("./data/real_lfc_data.csv")
cosimu::infer_dist(lfc_vector = real_lfc_data$log2FoldChange,save_yaml = "./data/real_asymmetrical_dist.yml")
cosimu::infer_dist(lfc_vector = real_lfc_data$log2FoldChange,save_yaml = "./data/real_symmetrical_dist.yml",sym_mode = TRUE)
