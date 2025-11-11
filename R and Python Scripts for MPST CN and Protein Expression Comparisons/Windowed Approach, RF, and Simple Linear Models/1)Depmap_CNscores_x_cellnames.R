# Load required packages
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ---------------------------
# Step 1: Read the CSV without headers
# ---------------------------
cn_raw <- read_csv(
  "data/raw/PortalOmicsCNGeneLog2.csv",
  col_names = FALSE,
  show_col_types = FALSE
)

# ---------------------------
# Step 2: Combine the first 9 rows to get the gene names vector
# ---------------------------
gene_rows <- cn_raw[1:9, -1]  # remove first column (likely empty)
gene_vector <- as.vector(t(gene_rows))  # flatten row-wise

# Clean gene names: remove leading commas, trim spaces, remove parentheses/HGNC IDs
gene_vector <- gene_vector %>%
  str_remove("^,") %>%
  str_trim() %>%
  str_remove("\\s*\\(.*\\)$")

# ---------------------------
# Step 3: Identify indices of genes of interest
# ---------------------------
genes_of_interest <- c("GPX4", "AIFM2", "MPST")
gene_indices <- which(gene_vector %in% genes_of_interest)

# ---------------------------
# Step 4: Extract numeric data (rows 10 onward)
# ---------------------------
data_rows <- cn_raw[-(1:9), ]  # rows 10 onward
numeric_data <- data_rows %>%
  select(1, all_of(gene_indices + 1))  # +1 because first column is ModelID

# Assign column names: first = ModelID, rest = genes of interest
colnames(numeric_data) <- c("CellID", gene_vector[gene_indices])

# ---------------------------
# Step 5: Pivot to long format
# ---------------------------
cn_long <- numeric_data %>%
  pivot_longer(
    cols = -CellID,
    names_to = "Gene",
    values_to = "Gene_CN_log2_value"
  )

# ---------------------------
# Step 6: Load Model metadata and merge
# ---------------------------
model_info <- read_csv(
  "data/raw/Model.csv",
  show_col_types = FALSE
)

cn_with_names <- cn_long %>%
  mutate(CellID = str_trim(CellID)) %>%
  left_join(model_info %>% select(ModelID, CellLineName), by = c("CellID" = "ModelID")) %>%
  select(CellID, CellLineName, Gene, Gene_CN_log2_value)

# ---------------------------
# Step 7: Save final CSV
# ---------------------------
write_csv(cn_with_names, "data/processed/Gene_CN_log2_with_names.csv")

# ??? Done! CSV contains: CellID | CellLineName | Gene | Gene_CN_log2_value

