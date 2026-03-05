# Figure 1B
# Testis specificity among different gene category
library(ggplot2)

# load data and calculate log2
Gene_expression<-read.csv("data/Dpse_RNAseq_FPKM.csv")
Gene_expression$Testis2carcass<-as.numeric(Gene_expression$Testis2carcass)

Gene_expression$log2TC<-log(Gene_expression$Testis2carcass,base = 2)
Dpse_gene_list<-read.csv("data/gene_list.csv")

# annotate gene category
Gene_expression$gene_type<-NA
Gene_expression$gene_type[Gene_expression$Gene %in% Dpse_gene_list$Dpse_gene[Dpse_gene_list$Type=="New duplicated gene"]] <-"New Dup."
Gene_expression$gene_type[Gene_expression$Gene %in% Dpse_gene_list$Dpse_gene[Dpse_gene_list$Type=="New duplicated parental gene"]] <-"Parent"
Gene_expression$gene_type[Gene_expression$Gene %in% Dpse_gene_list$Dpse_gene[Dpse_gene_list$Type=="Single copy gene"]] <-"Single copy"
Gene_expression$gene_type[Gene_expression$Gene %in% Dpse_gene_list$Dpse_gene[Dpse_gene_list$Type=="De novo gene"]] <-"De novo"

Gene_expression<-Gene_expression[!is.na(Gene_expression$gene_type),]

Gene_expression$gene_type<-factor(Gene_expression$gene_type,levels = c("Single copy","Parent","New Dup.","De novo"),ordered = T)

# plot
pdf("testis_specificity.pdf",height = 2.5,width = 1.5)
ggplot(Gene_expression,mapping=aes(x=gene_type,y=log2TC,colour=gene_type,fill=gene_type))+ 
  geom_violin(trim=FALSE,width=0.7)+
  geom_boxplot(outlier.colour = NA,width=0.2,position = position_dodge(0.9))+ 
  theme(legend.position="right")+
  theme_bw() +
  xlab("")+ylab("log2 T/C")+
  theme(legend.title=element_blank(),legend.position="right")+
  scale_y_continuous(limits = c(-9,18),breaks = c(-5,0,5,10,15))+
  scale_color_manual(values =  c("white","white","white","white"),guide = "none")+
  scale_fill_manual(values = c("#C63652","#FFB379","#FFB379","#3482B2"),guide="none")
dev.off()