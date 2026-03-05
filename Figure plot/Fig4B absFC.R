# Figure 4B
# Violin plot of foldchange for each gene

library(Seurat)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)

# load data
Dpse_gene <- read.csv("data/gene_list.csv")
Dpse.integrated <- readRDS(file = "data/Dpse.integrated.rds")

# calculate foldchange for all genes
comparisons <- list(
  C1 = list("Eu Spermatocytes-E1",  "Para Spermatocytes-E1"),
  C2 = list("Eu Spermatocytes-E2",  c("Para1 Spermatocytes-E2","Para2 Spermatocytes-E2")),
  C3 = list("Eu Spermatocytes-L",   c("Para1 Spermatocytes-L","Para2 Spermatocytes-L")),
  S1 = list("Eu Spermatids-E",      c("Para1 Spermatids-E","Para2 Spermatids-E")),
  S2 = list("Eu Spermatids-L",      "Para Spermatids-L")
)

FC_list <- imap(comparisons, function(comp, name) {
  
  res <- FindMarkers(
    Dpse.integrated,
    ident.1 = comp[[1]],
    ident.2 = comp[[2]],
    logfc.threshold = 0
  )
  
  tibble(
    Dpse_name = rownames(res),
    !!name := res$avg_log2FC
  )
})


FClist<-list(C1_FC,C2_FC,C3_FC,S1_FC,S2_FC)
FC_col<-c("C1","C2","C3","S1","S2")
for (i in 1:5) {
  FClist[[i]]<-cbind(rownames(FClist[[i]]),FClist[[i]])
  FClist[[i]]<-FClist[[i]][,c(1,3)]
  colnames(FClist[[i]])<-c("Dpse_name",FC_col[i])
}

Total_FC <- reduce(FC_list, full_join, by = "Dpse_name") %>%
  mutate(across(-Dpse_name, ~replace_na(.x, 0)))

# calculate maxFC
fc_cols <- setdiff(colnames(Total_FC), "Dpse_name")

Total_FC <- Total_FC %>%
  mutate(
    sumFC = rowSums(across(all_of(fc_cols))),
    maxFC = ifelse(
      sumFC >= 0,
      pmap_dbl(select(., all_of(fc_cols)), max),
      pmap_dbl(select(., all_of(fc_cols)), min)
    ),
    absFC = abs(maxFC)
  )


# classify genes
Total_FC <- Total_FC %>%
  mutate(
    type = case_when(
      Dpse_name %in% Dpse_gene$Dpse_gene[Dpse_gene$Type=="New duplicated gene"] ~ "New dup.",
      Dpse_name %in% Dpse_gene$Dpse_gene[Dpse_gene$Type=="New duplicated parental gene"] ~ "Parent",
      Dpse_name %in% Dpse_gene$Dpse_gene[Dpse_gene$Type=="De novo gene"] ~ "De novo",
      Dpse_name %in% Dpse_gene$Dpse_gene[Dpse_gene$Type=="Single copy gene"] ~ "Single copy",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(type))

Total_FC$type <- factor(
  Total_FC$type,
  levels = c("Single copy","Parent","New dup.","De novo"),
  ordered = TRUE
)


pdf("Vln_absFC.pdf",height = 2.5,width = 3)
ggplot(Total_FC, aes(x = type, y = absFC, fill = type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  theme_bw() +
  labs(x = NULL, y = "abs log2FC") +
  scale_y_continuous(limits = c(0,2), breaks = c(0,0.5,1,1.5,2)) +
  scale_fill_manual(
    values = c(
      "Single copy" = "#C63652",
      "Parent"      = "#FFB379",
      "New dup."    = "#FFB379",
      "De novo"     = "#3482B2"
    )
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )


