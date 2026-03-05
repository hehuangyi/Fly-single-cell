# Figure 3H
# Gradient gene heatmap

library(Seurat)
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)

# load data
Dpse.integrated <- readRDS(file = "data/Dpse.integrated.rds")
gradient_gene <- read.csv(file = "data/gradient.csv")

# calculate average expression
avg_exp <- AverageExpression(
  Dpse.integrated,
  assays = "RNA",
  slot = "data",
  group.by = "cell_types"
)

exp_mat <- as.data.frame(avg_exp$RNA)

cell_order <- c(
  "Eu Spermatocytes-E2","Para2 Spermatocytes-E2","Para1 Spermatocytes-E2",
  "Eu Spermatocytes-L","Para2 Spermatocytes-L","Para1 Spermatocytes-L",
  "Eu Spermatids-E","Para2 Spermatids-E","Para1 Spermatids-E"
)

exp_mat <- exp_mat[, cell_order]

# extract gradient genes
exp_mat <- exp_mat[rownames(exp_mat) %in% gradient_gene$Dpse_gene, ]

exp_mat <- exp_mat[match(gradient_gene$Dpse_gene, rownames(exp_mat)), ]

# z-score scaling
Dpse_expdata<-as.matrix(Dpse_expdata)
Dpse_expdata_scale<-t(scale(t(Dpse_expdata)))

# label eu/ para gradient
eu_genes   <- gradient_gene %>% filter(type == "Eu gradient") %>% pull(Dpse_gene)
para_genes <- gradient_gene %>% filter(type == "Para gradient") %>% pull(Dpse_gene)

exp_eu   <- exp_scale[rownames(exp_scale) %in% eu_genes, ]
exp_para <- exp_scale[rownames(exp_scale) %in% para_genes, ]

exp_scale_final <- rbind(exp_eu, exp_para)

row_split_factor <- factor(
  c(rep("Eu", nrow(exp_eu)),
    rep("Para", nrow(exp_para))),
  levels = c("Eu","Para")
)


# row annotation
mark_genes <- c(
  "rcd7","sowi","S-Lap1","loopin-1","S-Lap7","boly",
  "CG5048","Mst98Ca","Dpy-30L2","chic","nes",
  "CG16719","CP110"
)

mark_index <- which(rownames(exp_scale_final) %in% mark_genes)

ha <- rowAnnotation(
  link = anno_mark(
    at = mark_index,
    labels = rownames(exp_scale_final)[mark_index]
  )
)

# plot

pdf(output_file, height = 8, width = 4)

Heatmap(
  exp_scale_final,
  name = "Z-score",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2d77ab", "#f0ead2", "#e41a1c")
  ),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  column_split = rep(1:3, each = 3),
  row_split = row_split_factor,
  left_annotation = rowAnnotation(
    block = anno_block(
      gp = gpar(fill = c("#e41a1c", "#2d77ab"))
    )
  ),
  right_annotation = ha
)

dev.off()