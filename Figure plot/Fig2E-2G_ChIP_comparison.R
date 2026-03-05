# Figure 2E-2G
# IP/IN comparison among chromosomes

library(ggplot2)

# load files
all_chip_data <- read.csv("data/Gene_IPIN_value.csv")
gene_chr_file <- read.csv("data/Dpse_gene_chromosome.csv")

# define chromosome gene sets
AX_gene <- Dpse_gene$Gene[Dpse_gene$Chromosome=="AX"]
DX_gene <- Dpse_gene$Gene[Dpse_gene$Chromosome=="DX"]
Autosome_gene <- Dpse_gene$Gene[Dpse_gene$Chromosome %in% c("2","3","4")]

# markers
head_marker <- c("Hm_H3K9me2","Hm_H4K20me3","Hm_H4K16AC")
Testis_marker <-c ("Tm_H3K9me2","Tm_H4K20me3","Tm_H4K16AC")
ChIP_marker <-c ("H3K9me2","H4K20me3","H4K16AC")

# plotting function

plot_chip <- function(marker_index){
  
  chip_value <- all_chip_data[
    all_chip_data$Sample %in% c(head_marker[marker_index],
                                testis_marker[marker_index]), ]
  
  chip_value$chr <- NA
  chip_value$chr[chip_value$Gene %in% Autosome_gene] <- "Autosome"
  chip_value$chr[chip_value$Gene %in% AX_gene] <- "AX"
  chip_value$chr[chip_value$Gene %in% DX_gene] <- "DX"
  
  chip_value <- chip_value[!is.na(chip_value$chr), ]

  chip_value$chr <- factor(chip_value$chr,
                           levels = c("Autosome","AX","DX"))
  chip_value$Sample <- factor(chip_value$Sample,
                              levels = c(head_marker[marker_index],
                                         testis_marker[marker_index]))
  
  out_pdf <- paste0("Chip-seq/ChIP_comparison_",
                    chip_marker[marker_index], ".pdf")
  
  pdf(out_pdf, width = 4, height = 2)
  
  p <- ggplot(chip_value,
              aes(x = Sample, y = IP2IN, fill = chr, color = chr)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(outlier.colour = NA,
                 width = 0.2,
                 position = position_dodge(0.9)) +
    theme_bw() +
    ylab("IP/IN") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.title = element_blank()) +
    scale_y_continuous(limits = c(0,3), breaks = c(1,2,3)) +
    scale_x_discrete(labels = c("Head","Testis")) +
    scale_color_manual(values = c("white","white","white")) +
    scale_fill_manual(values = c("#3288bd","#C63652","#FFB379"))
  
  print(p)
  dev.off()
}

# run
for (i in seq_along(chip_marker)) {
  plot_chip(i)
}
