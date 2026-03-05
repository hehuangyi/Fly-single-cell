# Figure 2C
# Dpse AX-lined gene/Autosome expression ratio

library(Seurat)
library(ggplot2)
library(dplyr)

# load files
Dpse.integrated <- readRDS(file = "data/Dpse.integrated.rds")
gene_chr_file <- read.csv("data/Dpse_gene_chromosome.csv")
gene_list <- read.csv("data/gene_list.csv")

# define chromosome gene sets
AX_gene <- gene_chr_file$Gene[gene_chr_file$Chromosome=="AX"]
singlecopy_gene <- intersect(AX_gene,gene_list$Dpse_gene[gene_list$Type=="Single copy gene"])
newdup_gene <- intersect(AX_gene,gene_list$Dpse_gene[gene_list$Type=="New duplicated gene"])
denovo_gene<-intersect(AX_gene,gene_list$Dpse_gene[gene_list$Type=="De novo gene"])

# extract count matrix
count_mat <- Dpse.integrated@assays$RNA@counts

singlecopy_count <- colSums(count_mat[rownames(count_mat) %in% singlecopy_gene, ])
singlecopy_num <- sum(rownames(count_mat) %in% singlecopy_gene)

newdup_count <- colSums(count_mat[rownames(count_mat) %in% newdup_gene, ])
newdup_num <- sum(rownames(count_mat) %in% newdup_gene)

denovo_count <- colSums(count_mat[rownames(count_mat) %in% denovo_gene, ])
denovo_num <- sum(rownames(count_mat) %in% denovo_gene)

# Calculate normalized ratios
singlecopy_to_auto <- singlecopy_count / auto_counts / singlecopy_num * auto_num
newdup_to_auto <- newdup_counts / auto_counts / newdup_num * auto_num
denovo_to_auto  <- denovo_counts  / auto_counts / denovo_num  * auto_num
df_singlecopy <- data.frame(value = singlecopy_to_auto,
                    celltype = Dpse.integrated$cell_types,
                    stages = Dpse.integrated$stage,
                    gene = "Single copy")

df_newdup <- data.frame(value = newdup_to_auto,
                    celltype = Dpse.integrated$cell_types,
                    stages = Dpse.integrated$stage,
                    gene = "New Dup.")

df_denovo <- data.frame(value = denovo_to_auto,
                   celltype = Dpse.integrated$cell_types,
                   stages = Dpse.integrated$stage,
                   gene = "De novo")

sex_to_auto <- rbind(df_singlecopy, df_newdup, df_denovo)

# order factor
cell_order<- c("Spermatogonia-E","Spermatogonia-L",
               "Eu Spermatocytes-E1","Eu Spermatocytes-E2","Eu Spermatocytes-L","Eu Spermatids-E","Eu Spermatids-L",
               "Cyst-E","Cyst-L")
cell_type<- c("Spermatogonia-E","Spermatogonia-L",
              "Spermatocytes-E1","Spermatocytes-E2","Spermatocytes-L","Spermatids-E","Spermatids-L",
              "Cyst-E","Cyst-L")


sex_to_auto$celltype <- factor(sex_to_auto$celltype,
                               levels = cell_order,
                               ordered = TRUE)
sex_to_auto<-sex_to_auto[!is.na(sex_to_auto$celltype),]
sex_to_auto$gene <- factor(sex_to_auto$gene,
                           levels = c("Single copy","New Dup.","De novo")
)

# for Figure 2C, plot adult stage
sex_to_auto <- subset(sex_to_auto, stages == "Dpse adult")

# summary statistics
summary_df <- sex_to_auto %>%
  group_by(celltype, stages, gene) %>%
  summarise(
    median = median(value),
    sd = sd(value),
    .groups = "drop"
  )

#plot
pdf("results/figures/Dpse_AXgenetoAutosome_ratio.pdf", width = 6, height = 4)
ggplot(summary_df,
       aes(x = celltype, y = median,
           ymin = median - sd,
           ymax = median + sd,
           color = gene, fill = gene,
           group = gene)) +
  theme_bw() +
  geom_pointrange(size = 0.8, shape = 21) +
  geom_line() +
  scale_x_discrete(breaks = cell_order,
                   labels = cell_type) +
  scale_y_continuous(limits = c(0,0.9),
                     breaks = seq(0,0.9,0.1)) +
  scale_color_manual(values = c("#C63652","#FFB379","#3482B2")) +
  scale_fill_manual(values = c("#C63652","#FFB379","#3482B2")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  ylab("AX/Auto ratio") +
  xlab("")

dev.off()
