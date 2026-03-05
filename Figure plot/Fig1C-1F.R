# Figure 1C-1F

library(Seurat)
library(ggplot2)

Dpse.integrated<-readRDS("data/Dpse.integrated.rds")

# Plot Figure 1C
# Dot plot for cell type markers in D. pseudoobscura

Dpse_gene<-c("LOC4803214","LOC117183186","LOC4817921",
             "LOC26534190","LOC6900097","LOC6896840",
             "LOC4813925","LOC4803852",
             "LOC4817866","LOC4814025","LOC4801948",
             "LOC4802953","LOC6897136","LOC4816276",
             "LOC4818017","LOC4816969","LOC4801502")
Dmel_name<-c("bam","aub","vas","AGO3","aly","can","sa","CycB","twe",
             "bol","soti","wa-cup","f-cup","tj","eya","Fas3","MtnA")

pdf(file="cell_marker.pdf",width = 6.5,height = 4.5)

DotPlot(Dpse.integrated,features = Dpse_gene, 
        dot.scale = 3,scale.by = 'size' ,group.by = "cell_types") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1))+xlab("")+ylab("")+
  scale_colour_gradient2(low = "#2d77ab", mid="#f0ead2", high = "#e41a1c")+
  scale_x_discrete(breaks=Dpse_gene,labels=Dmel_name)

dev.off()

# Plot Figure 1D
# UMAP projection of D.pse cells
Dpse.integrated$cell_types<-
  factor(Dpse.integrated$cell_types,
         levels = c("Spermatogonia-E","Spermatogonia-L",
                    "Eu Spermatocytes-E1","Eu Spermatocytes-E2","Eu Spermatocytes-L","Eu Spermatids-E","Eu Spermatids-L",
                    "Para Spermatocytes-E1","Para1 Spermatocytes-E2","Para2 Spermatocytes-E2", "Para1 Spermatocytes-L","Para2 Spermatocytes-L","Para1 Spermatids-E","Para2 Spermatids-E","Para Spermatids-L",
                    "Cyst-E","Cyst-L","Hub cells","Epithelial cells"),ordered = T)


umapcolor<-c("#ffffbf","#fee08b",
             "#FFCD82","#FFB379","#FF6D68","#C63652","#9e0142",
             "#E5F0B0", "#D2F995","#66c2a5",
             "#4A9451","#7FA2DC","#3288bd","#5e4fa2","#115881",
             "#6A5A53","#BE9865","#9E899B","grey")

pdf("Dpse_umap.pdf",width = 8,height = 4)
DimPlot(Dpse.integrated,group.by = "cell_types",split.by = "stage",ncol = 2)+
  scale_color_manual(values=umapcolor)+labs(title = "")+
  coord_fixed()+NoAxes()
dev.off()

# Plot Figure 1E
# UMAP projection of integrated fly testis cells
fly.integrated<-readRDS("data/fly_testis_integrated.rds")

pdf("fly_testis_umap.pdf",width = 6,height = 4)
DimPlot(fly.integrated,group.by = "species",shuffle =T)+coord_fixed()+scale_color_manual(values =  c("#feb095","#8da8bc"))+NoAxes()+labs(title = "")
dev.off()

# Plot Figure 1F
# Cell type frequency plot

Dpse_freq<-data.frame(table(Dpse.integrated$cell_types,Dpse.integrated$stage))

Dpse_freq$Var1<-factor(Dpse_freq$Var1, levels =  c("Spermatogonia-E","Spermatogonia-L",
                                                   "Eu Spermatocytes-E1","Eu Spermatocytes-E2","Eu Spermatocytes-L","Eu Spermatids-E","Eu Spermatids-L",
                                                   "Para Spermatocytes-E1","Para1 Spermatocytes-E2","Para2 Spermatocytes-E2","Para1 Spermatocytes-L","Para2 Spermatocytes-L", "Para1 Spermatids-E","Para2 Spermatids-E","Para Spermatids-L",
                                                   "Cyst-E","Cyst-L","Hub cells","Epithelial cells"),ordered = T)

Dpse_freq<-Dpse_freq[!is.na(Dpse_freq$Var1),]

pdf("Dpse_cell_frequency.pdf",height = 5,width = 4.5)
ggplot(Dpse_freq, aes(x=Var2, y=Freq, group=Var1)) +
  geom_bar(stat = "identity", aes(fill=Var1),position = "fill")+ 
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 45, hjust=0.5, vjust=0.5), 
        legend.position= "right",legend.title = element_blank()) + 
  scale_fill_manual(values =umapcolor)+
  xlab("") +ylab("Cell number percentage")
dev.off()

