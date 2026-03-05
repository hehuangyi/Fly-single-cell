# Figure 4E
# Pie chart for functional categories

library(ggplot2)
library(dplyr)
library(tidyr)

# define data
type <- c("Testis silenced",
          "Neofunctionalization",
          "Cell-type conserved")

new_dup_nums     <- c(13, 60, 22)
single_copy_nums <- c(971, 2966, 6380)

df <- tibble(
  type = type,
  `New duplicate` = new_dup_nums,
  `Single copy`   = single_copy_nums
)

df_long <- df %>%
  pivot_longer(
    cols = -type,
    names_to = "Category",
    values_to = "Count"
  ) %>%
  group_by(Category) %>%
  mutate(
    Percentage = round(Count / sum(Count) * 100, 1),
    label = paste0(type, " (", Percentage, "%)")
  ) %>%
  ungroup()

# plot
my_colors <- c(
  "Testis silenced"      = "#3482B2",
  "Neofunctionalization" = "#C63652",
  "Cell-type conserved"  = "#FFB379"
)

plot_pie <- function(data, category_name) {
  
  subdata <- data %>% filter(Category == category_name)
  
  ggplot(subdata, aes(x = "", y = Count, fill = type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = my_colors,
      labels = subdata$label
    ) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 9)
    ) +
    ggtitle(category_name)
}

categories <- unique(df_long$Category)

for (cat in categories) {
  
  output_file <- file.path(
    output_dir,
    paste0("pie_function_", gsub(" ", "_", cat), ".pdf")
  )
  
  ggsave(
    output_file,
    plot = plot_pie(df_long, cat),
    width = 4,
    height = 2.5
  )
  
  message("Saved: ", output_file)
}
