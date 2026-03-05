# Figure 2A
# Dpse XY/Autosome expression ratio

library(Seurat)
library(ggplot2)
library(dplyr)

# load files
Dpse.integrated <- readRDS(file = "data/Dpse.integrated.rds")
gene_chr_file <- read.csv("data/Dpse_gene_chromosome.csv")

# define chromosome gene sets
AX_gene <- gene_chr_file$Gene[gene_chr_file$Chromosome=="AX"]
DX_gene <- gene_chr_file$Gene[gene_chr_file$Chromosome=="DX"]
Autosome_gene <- gene_chr_file$Gene[gene_chr_file$Chromosome %in% c("2","3","4")]
Y_gene <- gene_chr_file$Gene[gene_chr_file$Chromosome=="Y"]

# extract count matrix
count_mat <- Dpse.integrated@assays$RNA@counts

auto_counts <- colSums(count_mat[rownames(count_mat) %in% Autosome_gene, ])
auto_num <- sum(rownames(count_mat) %in% Autosome_gene)

AX_counts <- colSums(count_mat[rownames(count_mat) %in% AX_gene, ])
AX_num <- sum(rownames(count_mat) %in% AX_gene)

DX_counts <- colSums(count_mat[rownames(count_mat) %in% DX_gene, ])
DX_num <- sum(rownames(count_mat) %in% DX_gene)

y_counts <- colSums(count_mat[rownames(count_mat) %in% Y_gene, ])
y_num <- sum(rownames(count_mat) %in% Y_gene)

# Calculate normalized ratios
AX_to_auto <- AX_counts / auto_counts / AX_num * auto_num
DX_to_auto <- DX_counts / auto_counts / DX_num * auto_num
y_to_auto  <- y_counts  / auto_counts / y_num  * auto_num
df_AX <- data.frame(value = AX_to_auto,
                    celltype = Dpse.integrated$cell_types,
                    stages = Dpse.integrated$stage,
                    gene = "AX")

df_DX <- data.frame(value = DX_to_auto,
                    celltype = Dpse.integrated$cell_types,
                    stages = Dpse.integrated$stage,
                    gene = "DX")

df_y <- data.frame(value = y_to_auto,
                   celltype = Dpse.integrated$cell_types,
                   stages = Dpse.integrated$stage,
                   gene = "YD")

sex_to_auto <- rbind(df_AX, df_DX, df_y)

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
                         levels = c("AX","DX","YD")
                         )
# for Figure 2A, plot adult stage
sex_to_auto <- subset(sex_to_auto, stages == "D.pse adult")

# summary statistics
summary_df <- sex_to_auto %>%
  group_by(celltype, stages, gene) %>%
  summarise(
    median = median(value),
    sd = sd(value),
    .groups = "drop"
  )

#plot
pdf("results/figures/Dpse_XYtoAutosome_ratio.pdf", width = 6, height = 4)
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
  ylab("XY/Auto ratio") +
  xlab("")

dev.off()
