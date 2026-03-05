# Figure 4C
# Dot plot of gene proportion across stages

library(ggplot2)
library(dplyr)
library(tidyr)

# input
stage<-c("G1-G2","G2-C2","C1-C3","C3-S1","S2")
df <- tibble(
  stage = stage,
  `single copy` = c(49.3,18.3,10.6,23.2,4),
  parental      = c(33.3,20,16,29.3,1.3),
  child         = c(17.6,8.1,23,43.2,8.1),
  denovo        = c(14,8.5,10.9,55,11.6)
)

df_long <- df %>%
  pivot_longer(
    cols = -stage,
    names_to = "type",
    values_to = "proportion"
  ) %>%
  mutate(
    type = factor(
      type,
      levels = c("denovo","child","parental","single copy"),
      ordered = TRUE
    )
  )

# plot
pdf("dot_gene_stage.pdf",width = 6,height = 4)
ggplot(df_long, aes(x = stage, y = type)) +
  geom_point(
    aes(size = proportion, fill = type),
    shape = 21,
    color = "black"
  ) +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "denovo"      = "#3288bd",
      "child"       = "#FFB379",
      "parental"    = "#FFB379",
      "single copy" = "#C63652"
    )
  ) +
  scale_size(
    range = c(2, 20),
    name = "Proportion (%)"
  ) +
  labs(x = "Stage", y = NULL) +
  theme(
    legend.title = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )