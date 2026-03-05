############################################################
# Identification of Dpse-specific duplications
#
# This script:
# 1. Adds Muller element annotations to orthogroups
# 2. Identifies candidate duplication events
# 3. Identifies candidate translocation events
############################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(data.table)
library(purrr)

# ============================================================
# 1. Input files
# ============================================================
setwd("/Volumes/Elements/hhy/project/03.testis single cell new/00.supplementary/Code/New gene/")
ORTHOGROUP_FILE <- "Orthogroups.tsv"
GENECOUNT_FILE  <- "Orthogroups.GeneCount.tsv"
MULLER_DIR      <- "gene_muller/"

species_list <- c("Dmel", "Dmir", "Dpse", "Dsub", "Dvir", "Dwil")

# ============================================================
# 2. Helper functions
# ============================================================

# remove missing flag
rmN <- function(x) str_replace_all(x, "N", "")

# count gene copies in a cell
count_gene <- function(x) str_split(x, ", ") %>% unlist() %>% length()

# number of Muller elements
strlen <- function(x) nchar(x)

# ============================================================
# 3. Load orthogroups and append Muller
# ============================================================

orthogroups <- read_tsv(ORTHOGROUP_FILE)

for (sp in species_list) {
  
  message("Processing ", sp)
  
  muller_df <- read_csv(
    file.path(MULLER_DIR, paste0(sp, "_gene_muller.txt")),
    col_names = c("muller", "gene_id")
  )
  
  long_df <- tibble(
    orthogroup = orthogroups$Orthogroup,
    genes = orthogroups[[sp]]
  ) %>%
    separate_rows(genes, sep = ", ") %>%
    left_join(muller_df, by = c("genes" = "gene_id")) %>%
    group_by(orthogroup) %>%
    summarise(muller_concat = paste(na.omit(muller), collapse = ""),
              .groups = "drop")
  
  orthogroups <- left_join(
    orthogroups, long_df,
    by = c("Orthogroup" = "orthogroup")
  )
  
  colnames(orthogroups)[ncol(orthogroups)] <- paste0(sp, "_muller")
}

orthogroups[is.na(orthogroups)] <- "N"
raw_ortho <- copy(orthogroups)
orthogroups <- as.data.table(orthogroups)


# ============================================================
# PART I – DUPLICATIONS
# ============================================================
# ============================================================
# Scenario I – no missing Muller
# ============================================================

case0 <- orthogroups[
  rowSums(orthogroups[, c(2,4,6,7)] == "N") == 0
]

for (sp in c("Dwil","Dmir","Dmel","Dvir","Dsub")) {
  col <- paste0(sp, "_muller")
  case0[, (col) := map_chr(get(col), rmN)]
}

case0 <- case0[
  (Dmel_muller == Dvir_muller) &
    (Dmel_muller == Dwil_muller) &
    (Dpse_muller != Dmel_muller)
]

case0[, Dpse_n := map_chr(Dpse, count_gene)]
case0[, Dmel_gn := map_chr(Dmel, count_gene)]

dup_case0 <- case0[Dmel_gn == 1 & Dpse_n > 1]


# ============================================================
# Scenario II – allow one missing
# ============================================================

case1 <- raw_ortho %>%
  filter(rowSums(raw_ortho[, c(2,4,6,7)] == "N") == 1) %>%
  filter(Dpse_muller != "N") %>%
  filter(Dmel_muller != "N")

case1 <- as.data.table(case1)

for (sp in c("Dwil","Dmir","Dmel","Dvir","Dsub")) {
  col <- paste0(sp, "_muller")
  case1[, (col) := map_chr(get(col), rmN)]
}

case1[, Dmel_gn := map_chr(Dmel, count_gene)]
case1 <- case1[Dmel_gn == 1]

case1 <- case1[
  (Dmel_muller == Dvir_muller) |
    (Dmel_muller == Dwil_muller)
]

case1[, Dpse_n := map_chr(Dpse_muller, strlen)]
dup_case1 <- case1[Dpse_n > 1 & ! Dmel=="N"]


# ============================================================
# Extra duplications from gene count
# ============================================================

ortho_count <- read.table(GENECOUNT_FILE, header = TRUE)

extra <- ortho_count[
  ortho_count$Dmel == 1 &
    ortho_count$Dvir == 1 &
    ortho_count$Dwil == 1 &
    ortho_count$Dpse > 1, ]

extra <- raw_ortho[raw_ortho$Orthogroup %in% extra$Orthogroup, ]

extra <- as.data.table(extra)

for (sp in c("Dwil","Dmir","Dmel","Dvir","Dsub")) {
  col <- paste0(sp, "_muller")
  extra[, (col) := map_chr(get(col), rmN)]
}

extra[, Dmel_gn := map_chr(Dmel, count_gene)]
extra <- extra[Dmel_gn == 1]

extra <- extra[
  (Dmel_muller == Dvir_muller) |
    (Dmel_muller == Dwil_muller)
]

extra[, Dpse_n := map_chr(Dpse_muller, strlen)]


# final duplication set
all_duplications <- rbind(dup_case0, dup_case1, extra,fill=T)
all_duplications<-unique(all_duplications)
write.csv(
  all_duplications,
  file = "Dpse_duplications.csv",
  row.names = FALSE
)

############################################################
# PART II – TRANSLOCATIONS
############################################################

# ============================================================
# 0 missing
# ============================================================

case0 <- raw_ortho %>% 
  filter(rowSums(raw_ortho[, c(2,4,6,7)] == "N") == 0)

for (sp in c("Dwil","Dmir","Dmel","Dvir","Dsub")) {
  col <- paste0(sp, "_muller")
  case0[, (col) := map_chr(get(col), rmN)]
}

case0 <- case0[
  (Dmel_muller == Dvir_muller) &
    (Dmel_muller == Dwil_muller) &
    (Dpse_muller != Dmel_muller)
]

case0[, Dpse_n := map_chr(Dpse, count_gene)]
case0[, Dmel_gn := map_chr(Dmel, count_gene)]

trans_case0 <- case0[Dmel_gn == 1 & Dpse_n ==1]

# ============================================================
# 1 missing
# ============================================================

case1 <- raw_ortho %>%
  filter(rowSums(raw_ortho[, c(2,4,6,7)] == "N") == 1) %>%
  filter(Dpse_muller != "N") %>%
  filter(Dmel_muller != "N")

case1 <- as.data.table(case1)

for (sp in c("Dwil","Dmir","Dmel","Dvir","Dsub")) {
  col <- paste0(sp, "_muller")
  case1[, (col) := map_chr(get(col), rmN)]
}

case1[, Dmel_gn := map_chr(Dmel, count_gene)]
case1 <- case1[Dmel_gn == 1]

case1 <- case1[
  (Dmel_muller == Dvir_muller) |
    (Dmel_muller == Dwil_muller)
]

case1[, Dpse_n := map_chr(Dpse_muller, strlen)]
trans_case1 <- case1[Dpse_n == 1 & ! Dpse_muller == Dmel_muller & ! Dmel_muller==""]


# final translocation set
all_trans <- rbind(trans_case0, trans_case1)

csv(
  all_trans,
  file = "Dpse_translocations.csv",
  row.names = FALSE
)