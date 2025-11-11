
library(dplyr)
library(tidyr)
library(readr)

# ---------------------------
# Step 1: Load input files
# ---------------------------
prot <- read_csv("data/processed/combined_protein_erastin_resistance.csv")
cn <- read_csv("data/processed/Gene_CN_log2_with_names.csv")

# ---------------------------
# Step 2: Clean up column names
# ---------------------------
# Rename Erastin column for easier use
prot <- prot %>%
  rename(
    GPX4_protex = GPX4,
    AIFM2_protex = AIFM2,
    MPST_protex = MPST,
    Erastin_log2 = contains("ERASTIN"),
    CellLineName = `Cell Line Name`
  )

# ---------------------------
# Step 3: Pivot CN data wider
# ---------------------------
cn_wide <- cn %>%
  filter(Gene %in% c("MPST", "GPX4", "AIFM2")) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Gene_CN_log2_value,
    names_prefix = ""
  ) %>%
  rename(
    MPST_CN = MPST,
    GPX4_CN = GPX4,
    AIFM2_CN = AIFM2
  )

# ---------------------------
# Step 4: Merge by CellID
# ---------------------------
merged <- prot %>%
  left_join(cn_wide %>% select(CellID, CellLineName, MPST_CN, GPX4_CN, AIFM2_CN),
            by = "CellID")

# ---------------------------
# Step 5: Save final file
# ---------------------------
write_csv(merged, "data/processed/Merged_CN_Protein_complete.csv")

cat("??? Merged file saved to: data/processed/Merged_CN_Protein_complete.csv\n")
