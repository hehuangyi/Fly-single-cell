# Figure 1A
# Chromosome distribution among different gene sets

library(ggplot2)

# chromosome classes
chr_type<-c("Autosome","AX","DX")

# numbers for each gene class
singlecopy_nums<-c(6605,1622,2004)
newdup_nums<-c(78,53,12)
denovo_nums<-c(41,30,9)

# build dataframe
make_df <- function(nums, setname){
  df <- data.frame(chr_type = chr_type, nums = nums)
  df$set <- setname
  df$chr_type <- factor(df$chr_type, levels = chr_type)
  df$percent <- round(df$nums / sum(df$nums) * 100, 1)
  df$label <- paste0(df$chr_type, " (", df$percent, "%)")
  df
}

df_single <- make_df(singlecopy_nums, "Single-copy")
df_newdup <- make_df(newdup_nums, "New-dup")
df_denovo <- make_df(denovo_nums, "De novo")

# plot 
plot_pie <- function(df){
  ggplot(df, aes(x = "", y = nums, fill = chr_type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(title = unique(df$set), fill = "") +
    scale_fill_manual(values = c("#3482B2","#C63652","#FFB379"),
                      labels = df$label) +
    theme_void()
}

plot_pie(df_single)
plot_pie(df_newdup)
plot_pie(df_denovo)