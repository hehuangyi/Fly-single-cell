# generate D.mel and D.pse seurat object and process cross-species integration
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

set.seed(10)

# load ortholog list
Dpse_ortholog <- read.table("data/Dpse_ortholog.csv",
                            header = FALSE,
                            stringsAsFactors = FALSE)
ortholog_genes <- Dpse_ortholog$Dmel_gene

# creat D.mel Seurat objects
samples <- c("larvae1", "larvae2", "adult1", "adult2")
prefix  <- c("dmelL1-", "dmelL2-", "dmelA1-", "dmelA2-")

Dmel_list <- list()

for (i in seq_along(samples)) {
  
  message("Processing: ", samples[i])
  
  counts <- Read10X(data.dir = file.path("data/matrx/Dmel", samples[i]))
  
  colnames(counts) <- paste0(prefix[i], colnames(counts))
  
  # keep only ortholog genes
  counts <- counts[rownames(counts) %in% ortholog_genes, ]
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = paste0("Dmel_", samples[i]),
    min.cells = 3,
    min.features = 200
  )
  
  Dmel_list[[samples[i]]] <- obj

}

saveRDS(Dmel_list,
        file = file.path("results", "Dmel_list.rds"))


# load Dpse pre-processed data and remove doublet
samples <- c(
  "Dpse_adult1",
  "Dpse_adult2",
  "Dpse_larvae1",
  "Dpse_larvae2"
)

Dpse_list <- list()

for (sample in samples) {
  
  message("Processing: ", sample)
  
  obj <- readRDS(
    file.path("data/Dpse_doublet",
              paste0("doublet_", sample, ".rds"))
  )
  
  # Automatically find DF classification column
  df_col <- grep("DF.classifications",
                 colnames(obj@meta.data),
                 value = TRUE)
  
  obj <- subset(obj, obj[[df_col]] == "Singlet")
  
  obj$orig.ident <- sample
  
  Dpse_list[[sample]] <- obj
 
}

saveRDS(Dpse_list,
        file = file.path("results", "Dpse_list.rds"))

# map Dpse genes to ortholog table and rebuild Seurat objects
for (i in seq_along(Dpse.list)) {
  
  obj <- Dpse.list[[i]]
  obj <- subset(obj, nCount_RNA > 3000)
  
  counts <- GetAssayData(obj, slot = "counts")
  
  # match genes
  gene_index <- match(rownames(counts), Dpse_ortholog$Dpse_gene)
  valid_genes <- !is.na(gene_index)
  counts <- counts[valid_genes, ]
  
  # rename to ortholog gene name
  new_gene_names <- Dpse_ortholog$Dmel_gene[gene_index[valid_genes]]
  rownames(counts) <- new_gene_names
  
  # rebuild Seurat object
  Dpse_new_obj <- CreateSeuratObject(
    counts = counts,
    project = Dpse.ident[i]
  )
  
  Dpse_new_obj$orig.ident <- Dpse.ident[i]
  
  Dpse.mapped.list[[Dpse.ident[i]]] <- Dpse_new_obj
  
  )
}

saveRDS(
  Dpse.mapped.list,
  file = file.path("results", "Dpse_mapped_list.rds")
)


# integration D.mel and D.pse cells based on 1-1 orthologous genes
fly.list<-c(Dmel.list,Dpse.mapped.list)

for (i in 1:seq_along(fly.list)) {
  fly.list[[i]] <- NormalizeData(fly.list[[i]], verbose = FALSE)
  fly.list[[i]] <- FindVariableFeatures(fly.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

fly.anchors <- FindIntegrationAnchors(object.list = fly.list, dims = 1:30, k.anchor = 3)
fly.integrated <- IntegrateData(anchorset = fly.anchors, dims = 1:30)

# post-integration processing
DefaultAssay(fly.integrated) <- "integrated"
fly.integrated <- ScaleData(fly.integrated, verbose = FALSE)
fly.integrated <- RunPCA(fly.integrated, npcs = 30, verbose = FALSE)
fly.integrated <- RunUMAP(fly.integrated, reduction = "pca", dims = 1:30)
fly.integrated <- FindNeighbors(fly.integrated, reduction = "pca", dims = 1:30)
fly.integrated <- FindClusters(fly.integrated, resolution = 2)

saveRDS(fly.integrated,file = "data/fly_testis_integrated.rds")
