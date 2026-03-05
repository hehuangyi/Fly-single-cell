# Figure 3F
# GO bubble plot

library(ggplot2)
library(stringr)
library(dplyr)
library(readr)

# laod data
GO_filter <- read.csv("data//GO_summary_filter.csv")
GO_filter <- GO_filter %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(unique(Description)),
      ordered = TRUE
    ),
    
    # count gene number
    gene_num = str_count(Symbols, ",") + 1,
    
    # cap LogP at -10
    LogP = ifelse(LogP < -10, -10, LogP),
    
    # set type order
    Type = factor(
      Type,
      levels = c("DNB I", "DE I", "DNB II", "DE II", "Gradient"),
      ordered = TRUE
    )
  )

# add total gene numbers and calculate percentage 
total_gene_table <- tibble(
  Type = c("DE I", "DNB I", "DE II", "DNB II", "Gradient"),
  total_num = c(687, 487, 92, 395, 112)
)

GO_filter <- GO_filter %>%
  left_join(total_gene_table, by = "Type") %>%
  mutate(ratio = gene_num / total_num)

# plot
pdf(file = "GO_filter",width = 5.5,height = 3.5)
ggplot(GO_filter, aes(x = Type, y = Description)) +
  geom_point(aes(size = ratio, color = LogP)) +
  scale_colour_gradient2(
    low = "#C63652",
    mid = "#f0ead2",
    high = "#e41a1c"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Gene ratio",
    color = "LogP"
  )

dev.off()
