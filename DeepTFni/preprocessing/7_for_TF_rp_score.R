wd<-"/Your_path/"
library(Seurat)
Args <- commandArgs(trailingOnly = TRUE)
print(Args[1])
TF_file<-read.table(paste(wd,"Adjacency_matrix/",Args[1],".TF.list.txt",sep=''))
gene_rp_score<-Read10X_h5(paste(wd,Args[1],"/9_",Args[1],"_gene_score.h5",sep=''))
tf_rp_score = gene_rp_score[TF_file$V1,]
write.table(tf_rp_score,paste(wd,Args[1],"/10_",Args[1],"_TF_rp_score.txt",sep=''),sep="\t",quote=FALSE,row.names=FALSE)
