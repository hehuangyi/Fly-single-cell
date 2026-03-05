# Figure 4F
# Testis specificity of different functional categories

library(tidyverse)
library(ggplot2)

# load data
Dpse_stage<-read.csv("data/new_dup_functional_categories.csv")
Dpse_TC<-read.csv("data/Dpse_RNAseq_FPKM.csv")
Dmel_TC<-read.csv("data/Dmel_RNAseq_FPKM.csv")

Dpse_stage <- Dpse_stage %>%
  mutate(
    Functional.category = case_when(
      Functional.category == "Conserved cell types" ~ "Conserved",
      TRUE ~ Functional.category
    )
  )

# extract gene
extract_group <- function(category_name) {
  
  stage_sub <- Dpse_stage %>%
    filter(Functional.category == category_name)
  
  bind_rows(
    Dpse_TC %>%
      filter(Gene %in% stage_sub$New.Dup..gene) %>%
      mutate(source = "New dup"),
    
    Dpse_TC %>%
      filter(Gene %in% stage_sub$Parental.gene) %>%
      mutate(source = "Parental"),
    
    Dmel_TC %>%
      filter(Gene %in% stage_sub$Dmel.ortholog) %>%
      mutate(source = "Dmel")
    
  ) %>%
    mutate(type = category_name)
}

categories <- c(
  "Testis silenced",
  "Conserved",
  "Neofunctionalization"
)

all_ratio <- map_df(categories, extract_group)

# log2 transform
all_ratio <- all_ratio %>%
  filter(!is.na(Testis2carcass), Testis2carcass > 0) %>%   # avoid log2 error
  mutate(
    log2_ratio = log2(Testis2carcass),
    source = factor(source,
                    levels = c("New dup", "Parental", "Dmel"),
                    ordered = TRUE),
    type = factor(type,
                  levels = c("Testis silenced",
                             "Conserved",
                             "Neofunctionalization"),
                  ordered = TRUE)
  )

pdf("/Volumes/Elements/hhy/project/03.testis single cell new/16.stage/figure/vln_functionTC.pdf",width = 5,height = 3.5)
ggplot(all_ratio,
            aes(x = type,
                y = log2_ratio,
                fill = source,
                color = source)) +
  geom_violin(trim = FALSE, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2,
               outlier.shape = NA,
               position = position_dodge(0.9)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  ylab("log2(T/C)") +
  xlab("") +
  scale_y_continuous(limits = c(-7, 18)) +
  scale_color_manual(values = c("white","white","white"))+
  scale_fill_manual(values = c(
    "New dup"  = "#C63652",
    "Parental" = "#FFB379",
    "Dmel"     = "#3288bd"
  ))
dev.off()

