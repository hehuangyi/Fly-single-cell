# D.pse testis cell integration

library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)

set.seed(10)

# load list of Seurat object
Dpse.list<-readRDS("Dpse.list.rds")
Dpse.ident<-c("Dpse_adult1","Dpse_adult2","Dpse_larvae1","Dpse_larvae2")

# pre-filter cells
for (i in seq_along(Dpse.list)) {
  
  Dpse.list[[i]] <- subset(Dpse.list[[i]], nCount_RNA > 3000)
  Dpse.list[[i]]$orig.ident <- Dpse.list[[i]]@project.name
}

# normalize and find vairable features
for (i in seq_along(Dpse.list)) {
  Dpse.list[[i]] <- NormalizeData(Dpse.list[[i]], verbose = FALSE)
  Dpse.list[[i]] <- FindVariableFeatures(Dpse.list[[i]],
                                         selection.method = "vst",
                                         nfeatures = 2000,
                                         verbose = FALSE)
}

# integration
Dpse.anchors <- FindIntegrationAnchors(object.list = Dpse.list,
                                       dims = 1:30,
                                       k.anchor = 3)

Dpse.integrated <- IntegrateData(anchorset = Dpse.anchors, dims = 1:30)

# post-integration processing
DefaultAssay(Dpse.integrated) <- "integrated"
Dpse.integrated <- ScaleData(Dpse.integrated, verbose = FALSE)
Dpse.integrated <- RunPCA(Dpse.integrated, npcs = 30, verbose = FALSE)
Dpse.integrated <- RunUMAP(Dpse.integrated, reduction = "pca", dims = 1:30)
Dpse.integrated <- FindNeighbors(Dpse.integrated, reduction = "pca", dims = 1:30)
Dpse.integrated <- FindClusters(Dpse.integrated, resolution = 2)

# save
saveRDS(Dpse.integrated,file = "Dpse.integrated.rds")