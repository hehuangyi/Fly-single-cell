# Preprocessing and clustering for Dpse scRNA-seq
suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})


# load matrix
samples <- c(
  "Dpse_adult1",
  "Dpse_adult2",
  "Dpse_larvae1",
  "Dpse_larvae2"
)

data_list <- lapply(samples, function(x) {
  Read10X(data.dir = file.path("data", x, "outs/filtered_feature_bc_matrix"))
})

names(data_list) <- samples

# create seurat object
seurat_list <- list()

for (i in seq_along(samples)) {
  
  colnames(data_list[[i]]) <- paste0(samples[i], "-", colnames(data_list[[i]]))
  
  seurat_list[[i]] <- CreateSeuratObject(
    counts = data_list[[i]],
    project = samples[i],
    min.cells = 5,
    min.features = 500
  )
}


# QC and clustering

for (i in seq_along(seurat_list)) {
  obj <- seurat_list[[i]]

  count_q <- quantile(obj$nCount_RNA)
  count_iqr <- count_q[4] - count_q[2]
  
  feature_q <- quantile(obj$nFeature_RNA)
  feature_iqr <- feature_q[4] - feature_q[2]
  
  obj <- subset(
    obj,
    nCount_RNA < (count_q[4] + 1.5 * count_iqr) &
      nCount_RNA > (count_q[2] - 1.5 * count_iqr) &
      nFeature_RNA < (feature_q[4] + 1.5 * feature_iqr) &
      nFeature_RNA > (feature_q[2] - 1.5 * feature_iqr)
  )
  
  obj <- SCTransform(obj, verbose = TRUE, return.only.var.genes = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0)
  obj <- RunUMAP(obj, dims = 1:30)
  
  saveRDS(obj, file = file.path("results", paste0(samples[i], ".rds")))
  
  seurat_list[[i]] <- obj
}
