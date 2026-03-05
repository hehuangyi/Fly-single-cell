# Doublet detection using DoubletFinder

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(foreach)
  library(doParallel)
  library(dplyr)
})

samples <- c(
  "Dpse_adult1",
  "Dpse_adult2",
  "Dpse_larvae1",
  "Dpse_larvae2"
)

doublet_rate <- c(0.061, 0.061, 0.069, 0.061)

# Parallel setup
cl <- makeCluster(8)
registerDoParallel(cl)

foreach(i = seq_along(samples), .packages = c("Seurat","DoubletFinder","dplyr")) %dopar% {
  
  obj <- readRDS(file.path("results", paste0(samples[i], ".rds")))
  
  sweep.res <- paramSweep_v3(obj, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  nExp_poi <- round(doublet_rate[i] * ncol(obj))
  
  obj <- doubletFinder_v3(
    obj,
    PCs = 1:30,
    pN = 0.25,
    pK = mpK,
    nExp = nExp_poi,
    reuse.pANN = FALSE,
    sct = TRUE
  )
  
  saveRDS(obj, file.path("results", paste0("doublet_", samples[i], ".rds")))
}

stopCluster(cl)