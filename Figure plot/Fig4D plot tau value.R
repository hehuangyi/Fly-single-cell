# Figure 4D
# plot tau value for each gene expression stage

library(dplyr)
library(readr)
library(ggplot2)
library(tibble)

# load data
Dpse_stage<-read.csv("data/Dpse_stage.csv")
Dpse_fpkm <- read.csv("data/Dpse_FPKM.csv",row.names = 1)

# calculate tau, code is from https://github.com/severinEvo/gene_expression
tau_function <- function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t <- sum(1-x/max(x))/(length(x)-1)
}


Dpse_tau <- apply(as.matrix(Dpse_fpkm), 1, tau_function)

Dpse_tau <- data.frame(
  Dpse_gene = names(Dpse_tau),
  tau = Dpse_tau,
  row.names = NULL
)


# merge stage information
Dpse_tau <- Dpse_stage %>%
  left_join(Dpse_tau, by = "Dpse_gene") %>%
  filter(!is.na(tau))


# calculate median and sd
Dpse_tau_summary <- Dpse_tau %>%
  group_by(Stage) %>%
  summarise(
    median_tau = median(tau),
    sd_tau = sd(tau),
    .groups = "drop"
  ) %>%
  mutate(
    Stage = factor(
      Stage,
      levels = c("G1-G2","G2-C2","C1-C3","C3-S1","S2"),
      ordered = TRUE
    )
  )


pdf("tau_in_each_stage",width = 4,height = 3.2)
ggplot(Dpse_tau_summary, aes(x = Stage, y = median_tau)) +
  geom_pointrange(
    aes(
      ymin = median_tau - sd_tau,
      ymax = median_tau + sd_tau
    ),
    color = "#C63652",
    fill  = "#C63652",
    shape = 21,
    size = 1
  ) +
  geom_line(group = 1, color = "#C63652") +
  theme_bw() +
  labs(
    x = NULL,
    y = "Tau"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
dev.off()