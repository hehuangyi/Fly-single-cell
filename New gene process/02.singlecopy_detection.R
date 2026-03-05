############################################################
# Script: Extract Dpse single-copy orthologs
############################################################
# Load libraries
library(tidyverse)

# ============================================================
# 1. Input files
# ============================================================
ORTH_GROUP_FILE   <- "Orthogroups.tsv"
ORTH_COUNT_FILE   <- "Orthogroups.GeneCount.tsv"
OUTPUT_FILE       <- "Dpse_singlecopy.csv"

# ============================================================
# 2. Load data
# ============================================================
orth_group <- read_tsv(ORTH_GROUP_FILE)
orth_count <- read.table(ORTH_COUNT_FILE, header = TRUE)
dpse_convert <- read.table(DPSE_MAP_FILE, header = TRUE)

# ============================================================
# 3. Identify single-copy orthogroups
# ============================================================
single_orth <- orth_count$Orthogroup[
  orth_count$Dpse == 1 &
    orth_count$Dsub == 1 &
    orth_count$Dmel == 1
]

# Extract Dmel and Dpse gene IDs for these orthogroups
Dpse_singlecopy <- tibble(
  Dmel = orth_group$Dmel[orth_group$Orthogroup %in% single_orth],
  Dpse = orth_group$Dpse[orth_group$Orthogroup %in% single_orth]
)

# Remove species prefix from Dpse IDs
Dpse_singlecopy <- Dpse_singlecopy %>%
  mutate(Dpse = str_replace(Dpse, "Dpse_", ""))

# 4. Save output
# ============================================================
write.csv(Dpse_singlecopy, file = OUTPUT_FILE, row.names = FALSE)