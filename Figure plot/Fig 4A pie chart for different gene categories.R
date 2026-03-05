# Figure 4A
# Pie charts of gene categories

library(ggplot2)
library(dplyr)
library(tidyr)

# input data,the percentage of genes in each category
type <- c("Others","DEI","DEII","Gradient","DNBI","DNBII")

all_pct    <- c(86.9, 8.9, 1.6, 1.6, 6.2, 5.3)
all_dup    <- c(71.3,17.2, 6.6, 8.2,17.2,13.9)
all_denovo <- c(47.8,34.8,13,17.4,26.1,34.8)

df <- data.frame(
  type = type,
  All = all_pct,
  Dup = all_dup,
  De_novo = all_denovo
)

df_long <- df %>%
  pivot_longer(
    cols = c(All, Dup, De_novo),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    type = factor(type, levels = type, ordered = TRUE),
    label = paste0(type, " (", Percentage, "%)")
  )

my_colors <- c(
  "Others"   = "grey",
  "DEI"      = "#9e0142",
  "DEII"     = "#FF6D68",
  "Gradient" = "#FFB379",
  "DNBI"     = "#3288bd",
  "DNBII"    = "#66c2a5"
)

# plot
plot_pie <- function(data, category_name) {
  
  p <- ggplot(
    data = data %>% filter(Category == category_name),
    aes(x = "", y = Percentage, fill = type)
  ) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = my_colors,
      labels = unique(data$label[data$Category == category_name])
    ) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    ) +
    ggtitle(category_name)
  
  return(p)
}


categories <- c("All","Dup","De_novo")

for (cat in categories) {
  
  output_file <- file.path(
    "figure",
    paste0("pie_", cat, ".pdf")
  )
  
  ggsave(
    output_file,
    plot = plot_pie(df_long, cat),
    width = 4,
    height = 4
  )
  
  message("Saved: ", output_file)
}