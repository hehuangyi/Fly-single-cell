# Figure 2D
# Expression of new genes supporting MSCI

library(Seurat)
library(ggplot2)
library(cowplot)

Dmel.integrated <- readRDS("data/Dmel_integrated.rds")
Dmel_A <- subset(Dmel.integrated,stage=="D.mel adult")

Dpse.integrated<-readRDS("data/Dpse.integrated.rds")
Dpse_A<-subset(Dpse.integrated,stage=="D.pse adult")

# gene list
Dpse_dup_gene <- c("LOC6903594","LOC6897142","LOC6902854")
Dpse_parent_gene <- c("LOC4812761","LOC6900266","LOC4814556")
Dmel_gene <- c("polo","Rcd7","Arp2")
gene_order <- c("polo","Rcd7","Arp2")

gene_table <- data.frame(
  Dmel = rep(Dmel_gene, 2),
  Dpse = c(Dpse_dup_gene, Dpse_parent_gene)
)

# generate data for dot plot
Dmel.data <- DotPlot(Dmel_A,features = Dmel_gene,group.by = "cell_types")$data
Dpse.dup.data <- DotPlot(Dpse_A,features = Dpse_dup_gene,group.by = "cell_types")$data
Dpse.parent.data <- DotPlot(Dpse_A,features = Dpse_parent_gene,group.by = "cell_types")$data

Dmel.data$gene<-"D.mel"
Dpse.child.data$gene<-"New dup"
Dpse.pare.data$gene<-"Parental"

all.data <- rbind(Dpse.child.data,Dmel.data,Dpse.pare.data)
all.data <- merge(all.data, gene_table,by.x = "features.plot",by.y = "Dpse",all.x = T)
all.data$Dmel[is.na(all.data$Dmel)] <- as.character(all.data$features.plot[is.na(all.data$Dmel)])
all.data$gene <- factor(all.data$gene,levels =c("New dup","Parental","D.mel"),ordered = T)
all.data$Dmel <- factor(all.data$Dmel,levels = Dmel_gene,ordered = T)
all.data<-all.data[order(all.data$V1,all.data$gene),]

all.data$features.plot <- factor(all.data$features.plot,levels = unique(all.data$features.plot),ordered = T)
all.data <- all.data[all.data$id %in% cell_filter,]
all.data$id <- factor(all.data$id,
                      levels = rev(c("Spermatogonia-E","Spermatogonia-L",
                                     "Eu Spermatocytes-E1","Eu Spermatocytes-E2","Eu Spermatocytes-L","Eu Spermatids-E","Eu Spermatids-L")),ordered = T)
all.data$Dmel<-as.character(all.data$Dmel)

# set scale factor
all.data$scale_factor<-all.data$Dmel
all.data$scale_factor[all.data$gene %in% c("Parental","New dup")]<-paste0("Dpse_",all.data$scale_factor[all.data$gene %in% c("Parental","New dup")])

avg.exp.scaled <- sapply(
  unique(all.data$scale_factor),
  function(x){
    data.use <- all.data[all.data$scale_factor == x, "avg.exp"]
    data.use <- scale(data.use)
    data.use <- MinMax(data.use, min = -2.5, max = 2.5)
    return(data.use)
  }
)

all.data$avg.exp.scaled.gene <-  as.numeric(unlist(avg.exp.scaled))



ggplot(data = all.data, mapping = aes_string(x = 'features.plot', y = 'id')) +
  geom_point(mapping = aes_string(size = 'pct.exp',color = 'avg.exp.scaled.gene')) +
  guides(size = guide_legend(title = 'Percent Expressed')) +
  labs(x = '',y = '')  +
  scale_colour_gradient2(low = "#2d77ab", mid="#f0ead2", high = "#e41a1c",midpoint = 0.8)+
  theme_cowplot()+
  labs(x="",y="")


# bar plot for testis specificity
Dmel_data <- read.csv(file = "data/Dmel_RNAseq_FPKM.csv")
Dpse_data <- read.csv(file = "data/Dpse_RNAseq_FPKM.csv")

Dmel_TC <- Dmel_data[Dmel_tes_ratio$Gene %in% Dmel_gene,]
Dmel_TC$log2 <- log(Dmel_TC$Testis2carcass,base = 2)
Dmel_TC$species <- "D.mel"

Dpse_dup_TC <- Dpse_data[Dpse_data$Gene %in% Dpse_dup_gene,]
Dpse_dup_TC$log2 <- log(Dpse_child_ratio$Testis2carcass,base = 2)
Dpse_dup_TC$species <- "New dup"

Dpse_parent_TC <- Dpse_data[Dpse_data$Gene %in% Dpse_parent_gene,]
Dpse_parent_TC$log2 <- log(Dpse_parent_TC$Testis2carcass,base = 2)
Dpse_parent_TC$species <- "Parental"

all_TC <- rbind(Dmel_TC,Dpse_dup_TC,Dpse_parent_TC)
all_TC <- merge(all_TC, gene_table,by.x = "Gene",by.y = "Dpse",all.x = T)
all_TC$Dmel[is.na(all_TC$Dmel)] <- as.character(all_ratio$Gene[is.na(all_TC$Dmel)])
all_TC$species <- factor(all_TC$species,levels = c("New dup","Parental","D.mel"),ordered = T)
all_TC$Dmel < -factor(all_TC$Dmel,levels =Dmel_gene)
all_TC <- all_TC[order(all_TC$Dmel,all_TC$species),]

pdf("Bar_testis_specificity.pdf",width = 5,height = 2)
ggplot()+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.border = element_blank())+
  scale_fill_manual(values = c("#C63652","#3288bd","grey"))+ 
  geom_bar(data =all_ratio, aes(x = all_TC, y = log2,fill=species), stat = "identity")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_y_continuous(limits = c(-9,11),breaks = c(-5,0,5,10))+
  xlab("")+ylab("log2 T/C")
dev.off()


# bar plot for testis FPKM expression
Dmel_exp <- Dmel_data[Dmel_data$Gene %in% Dmel_gene,]
Dpse_dup_exp <- Dpse_data[Dpse_data$Gene %in% c(Dpse_dup_gene),]
Dpse_parent_exp <- Dpse_data[Dpse_data$Gene %in% c(Dpse_parent_gene),]

Dmel_exp$log2 <- log(Dmel_exp$M_testis_FPKM+1,base = 2)
Dmel_exp$species <- "D.mel"
Dpse_dup_exp$log2 <- log(Dpse_dup_exp$M_testis_FPKM+1,base = 2)
Dpse_dup_exp$species <- "New dup"
Dpse_parent_exp$log2 <- log(Dpse_parent_exp$M_testis_FPKM+1,base = 2)
Dpse_parent_exp$species <- "Parental"

all_exp <- rbind(Dmel_exp,Dpse_dup_exp,Dpse_parent_exp)
all_exp <- merge(all_exp, gene_table,by.x = "Gene",by.y = "Dpse",all.x = T)
all_exp$Dmel[is.na(all_exp$Dmel)] <- as.character(all_exp$Gene[is.na(all_exp$Dmel)])
all_exp$species <- factor(all_exp$species,levels = c("New dup","Parental","D.mel"),ordered = T)
all_exp$Dmel <- factor(all_exp$Dmel,levels = Dmel_gene,ordered = T)
all_exp <- all_exp[order(all_exp$Dmel,all_exp$species),]

pdf("Bar_testis_FPKM.pdf",width = 5,height = 2)
ggplot()+
  theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),panel.border = element_blank())+
  scale_fill_manual(values = c("#C63652","#3288bd","gray"))+ 
  geom_bar(data =all_exp, aes(x = X, y = log2,fill=species), stat = "identity")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_y_continuous(limits = c(0,10),breaks = c(0,5,10))+
  xlab("")+ylab("log2(FPKM+1)")
dev.off()
