# Figure 3D
# Expresion of meiotic arrest genes and mitochondiral complex genes

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)


# load data
Dmel.integrated <- readRDS("data/Dmel_integrated.rds")
Dpse.integrated <- readRDS("data/Dpse_integrated.rds")


# combine para1 and para2
Dpse.integrated$cell_types <- recode(
  Dpse.integrated$cell_types,
  "Para2 Spermatocytes-E2" = "Para1 Spermatocytes-E2",
  "Para2 Spermatocytes-L"  = "Para1 Spermatocytes-L",
  "Para2 Spermatids-E"     = "Para1 Spermatids-E"
)

Dpse_A <- subset(Dpse.integrated, stage == "D.pse Adult")
Dmel_A <- subset(Dmel.integrated, stage == "D.mel Adult")

# select genes to plot
Dmel_gene <- c("achi","CG5204","mip40","ATPsynCF6","mt:CoII")
Dpse_gene <- c("LOC4805467","LOC4816438","LOC4804612","LOC4800643","MT-COX2")

gene_table <- data.frame(
  Dmel_gene = Dmel_gene,
  Dpse_gene = Dpse_gene,
  stringsAsFactors = FALSE
)

# define cell types and cell order
cell_filter <- c("Spermatogonia-E","Spermatogonia-L",
               "Eu Spermatocytes-E1","Eu Spermatocytes-E2","Eu Spermatocytes-L","Eu Spermatids-E","Eu Spermatids-L",
               "Para Spermatocytes-E1","Para1 Spermatocytes-E2","Para1 Spermatocytes-L","Para1 Spermatids-E","Para Spermatids-L")

# extract dotplot data
extract_dotplot_data <- function(seurat_obj, features, species_label){
  
  dot <- DotPlot(
    seurat_obj,
    features = features,
    group.by = "cell_types"
  )
  
  df <- dot$data
  df$species <- species_label
  
  return(df)
}

Dmel.data <- extract_dotplot_data(Dmel_A, Dmel_gene, "D.mel")
Dpse.data <- extract_dotplot_data(Dpse_A, Dpse_gene, "D.pse")

# merge data
all.data <- bind_rows(Dpse.data, Dmel.data)
all.data <- left_join(
  all.data,
  gene_table,
  by = c("features.plot" = "Dpse_gene")
)

all.data$Dmel_gene[is.na(all.data$Dmel_gene)] <-
  all.data$features.plot[is.na(all.data$Dmel_gene)]

all.data <- all.data %>%
  filter(id %in% cell_filter) %>%
  mutate(
    species = factor(species,
                     levels = c("D.pse","D.mel"),
                     ordered = TRUE),
    id = factor(id,
                levels = rev(cell_filter),
                ordered = TRUE),
    Dmel_gene = factor(Dmel_gene,
                       levels = Dmel_gene,
                       ordered = TRUE),
    features.plot = factor(features.plot,
                           levels = unique(features.plot),
                           ordered = TRUE)
  ) %>%
  arrange(Dmel_gene, species)

# plot
pdf("meiotic_and_mt_gene_exp.pdf",width = 10,height = 4.5)
dot_plot <- ggplot(
  all.data,
  aes(x = features.plot, y = id)
) +
  geom_point(
    aes(size = pct.exp,
        color = avg.exp.scaled)
  ) +
  scale_colour_gradient2(
    low = "#2d77ab",
    mid = "#f0ead2",
    high = "#e41a1c",
    midpoint = 0.8
  ) +
  guides(size = guide_legend(title = "Percent Expressed")) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )
dev.off()