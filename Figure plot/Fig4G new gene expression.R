# Figure 4G
# Bubble plot for selected new gene expression

library(Seurat)
library(ggplot2)
library(cowplot)

# load data
Dmel.integrated<-readRDS("data/Dmel.integrated.rds")
Dpse.integrated<-readRDS("data/Dpse.integrated.rds")

Dpse$cell_types <- recode(
  Dpse$cell_types,
  "Para2 Spermatocytes-E3" = "Para1 Spermatocytes-E2",
  "Para2 Spermatocytes-E2" = "Para1 Spermatocytes-E2",
  "Para2 Spermatocytes-L"  = "Para1 Spermatocytes-L",
  "Para2 Spermatids-E"     = "Para1 Spermatids-E"
)


# gene list
Dmel_gene<-c(
  "RpS2","Rcd7","UQCR-C2","Arp2","fan","nes","CG13164","tbrd-1","eIF4E3","polo","Cp110","S-Lap5")
Dpse_C<-c(
  "LOC6901943","LOC6897142","LOC6897314","LOC6902854","LOC6902496","LOC6903581","LOC6901565","LOC6901771","LOC6902494","LOC6903594","LOC117183975","LOC6898081")
Dpse_P<-c(
  "LOC4816899","LOC6900266","LOC4812382","LOC4814556","LOC6900207","LOC4811983","LOC4803671","LOC6897274","LOC4814110","LOC4812761","LOC4815853","LOC4805168")

ene_table <- tibble(
  Dmel_gene = rep(Dmel_gene, 2),
  Dpse_gene = c(Dpse_C, Dpse_P),
  source    = rep(c("D.pse C", "D.pse P"), each = length(Dmel_gene))
)


# cell type and order
cell_filter<-c("Spermatogonia-E","Spermatogonia-L",
               "Eu Spermatocytes-E1","Eu Spermatocytes-E2","Eu Spermatocytes-L","Eu Spermatids-E","Eu Spermatids-L",
               "Para Spermatocytes-E1","Para1 Spermatocytes-E2","Para1 Spermatocytes-L","Para1 Spermatids-E","Para Spermatids-L")

# extract dotplot data 
ene_table <- tibble(
  Dmel_gene = rep(Dmel_gene, 2),
  Dpse_gene = c(Dpse_C, Dpse_P),
  source    = rep(c("D.pse C", "D.pse P"), each = length(Dmel_gene))
)

all_data <- bind_rows(Dmel_data, DpseC_data, DpseP_data)


all_data <- all_data %>%
  left_join(gene_table,
            by = c("features.plot" = "Dpse_gene")) %>%
  mutate(
    Dmel_gene = ifelse(is.na(Dmel_gene),
                       as.character(features.plot),
                       Dmel_gene)
  )

all_data <- all_data %>%
  filter(id %in% cell_filter) %>%
  mutate(
    id = factor(id, levels = rev(cell_filter)),
    species = factor(species,
                     levels = c("D.pse C","D.pse P","D.mel")),
    Dmel_gene = factor(Dmel_gene,
                       levels = Dmel_gene),
    features.plot = factor(features.plot,
                           levels = unique(features.plot))
  ) %>%
  arrange(Dmel_gene, species)

# plot
pdf("new_gene_exp.pdf",width = 14,height = 4.5)

ggplot(all_data,
       aes(x = features.plot,
           y = id)) +
  geom_point(aes(size = pct.exp,
                 color = avg.exp.scaled)) +
  scale_size(range = c(0.5, 3),
             limits = c(0, 100),
             name = "Percentage") +
  scale_colour_gradient2(
    low = "#2d77ab",
    mid = "#f0ead2",
    high = "#e41a1c",
    midpoint = 0.8
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )
dev.off()
