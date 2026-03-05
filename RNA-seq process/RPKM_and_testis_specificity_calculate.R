# Calculate FPKM and testis specificity score for RNA-seq

library(dplyr)
# load data and extract gene length
gene_counts<-read.csv("RNA-seq/Gene_RNAseq_counts.csv",row.names = 1,stringsAsFactors = FALSE)
gene_length<-as.numeric(gene_counts$Length)
gene_counts<-gene_counts[,-1]
gene_counts[] <- lapply(gene_counts, as.numeric)

# calculate FPKM 
kb <- gene_length/ 1000
rpk <- sweep(gene_counts, 1, kb, "/")
fpkm <- sweep(rpk, 2, colSums(gene_counts) / 1e6, "/")

# get average expression
gene_expression <- data.frame(
  M_carcass_FPKM = rowMeans(fpkm[, c(1, 2)]),
  M_testis_FPKM  = rowMeans(fpkm[, c(3, 4)])
)

# calculate testis specificity score
gene_expression$Testis2carcass<-(gene_expression$M_testis_FPKM+1e-6)/(gene_expression$M_carcass_FPKM+1e-6)

# output
write.csv(gene_expression,file = "RNA-seq/Gene_RNAseq_FPKM.csv")


