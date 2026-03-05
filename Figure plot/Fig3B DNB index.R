# Figure 3B
# DNB index

library(ggplot2)
library(ggpubr)

# load data
Dpse_DNB <- read.csv("data/Dpse_adult_DNBindex.csv")
Dmel_DNB <- read.csv("data/Dmel_adult_DNBindex.csv")

# define cell types
cell_type <- c("G1","G2","eu-C1","eu-C2","eu-C3","eu-S1","eu-S2",
             "para-C1","para-C21","para-C31","para-S11","para-C22","para-C32","para-S12","para-S2")

eu <- c("G1","G2","eu-C1","eu-C2","eu-C3","eu-S1","eu-S2")
para1 <- c("G1","G2","para-C1","para-C21","para-C31","para-S11","para-S2")
para2 <- c("G1","G2","para-C1","para-C22","para-C32","para-S12","para-S2")

# function: calculate mean and sd from even columns
calculate_DNB_stats <- function(df, cell_labels, cell_line_name){
  even_cols <- seq(2, ncol(df), by = 2)
  mean_vals <- colMeans(df[, even_cols])
  sd_vals   <- apply(df[, even_cols], 2, sd)
  
  data.frame(
    cell_type = cell_labels,
    DNB       = mean_vals,
    sd        = sd_vals,
    cell_line = cell_line_name
  )
}

# calculate Dpse stats
Dpse_stats <- calculate_DNB_stats(
  Dpse_DNB,
  cell_type,
  "Dpse_all"
)

Dpse_eu    <- subset(Dpse_stats, cell_type %in% eu)
Dpse_eu$cell_line <- "eu"
Dpse_eu$x <- 1:7

Dpse_para1 <- subset(Dpse_stats, cell_type %in% para1)
Dpse_para1$cell_line <- "para1"
Dpse_para1$x <- 1:7

Dpse_para2 <- subset(Dpse_stats, cell_type %in% para2)
Dpse_para2$cell_line <- "para2"
Dpse_para2$x <- 1:7

# calculate Dmel stats
Dmel_stats <- calculate_DNB_stats(
  Dmel_DNB,
  "Dmel",
  cell_type[1:7]
)
Dmel_stats$x <- 1:7
Dmel_stats$cell_line<-"Dmel"

# plot
DNB_allspecies <- rbind(
  Dpse_eu,
  Dpse_para1,
  Dpse_para2,
  Dmel_stats
)

DNB_allspecies$x <- factor(DNB_allspecies$x)

DNB_plot <- ggplot(
  DNB_allspecies,
  aes(x = x, y = DNB, group = cell_line)
) +
  theme_bw() +
  geom_line(aes(linetype = cell_line,
                color = cell_line)) +
  geom_pointrange(
    aes(ymin = DNB - sd,
        ymax = DNB + sd,
        color = cell_line),
    size = 0.8
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  scale_color_manual(
    values = c("grey","#cf3553","#3288bd","#5E4FA2")
  ) +
  scale_x_discrete(
    breaks = 1:7,
    labels = c("G1","G2","C1","C2","C3","S1","S2")
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 0.8, by = 0.1),
    limits = c(0.2, 0.85)
  ) +
  labs(
    x = "",
    y = "LDNB Index (Normalized)",
    title = ""
  )

pdf("DNB_index.pdf",width = 7 ,height = 3)

DNB_plot

dev.off()