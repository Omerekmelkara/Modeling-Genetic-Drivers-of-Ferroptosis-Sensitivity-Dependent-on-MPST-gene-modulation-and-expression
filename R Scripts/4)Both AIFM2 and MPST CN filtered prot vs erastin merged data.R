library(dplyr)
library(readr)

# === Step 1: Load merged proteomics + Erastin data ===
merged <- read_csv("data/processed/Merged_Proteomics_Erastin.csv")

# === Step 2: Load unique gene list from CN segments ===
unique_genes <- read_csv("data/processed/Unique_Genes.csv") %>%
  pull(Gene)   # vector of gene names

# build z-score column names
gene_cols <- paste0(unique_genes, "_z")

# keep only the ones that exist in merged data
gene_cols <- intersect(gene_cols, colnames(merged))

cat("??? Found", length(gene_cols), "genes from CN list present in merged dataset\n")

# === Step 3: Filter by AIFM2_z > 1 ===
filtered <- merged %>%
  filter(AIFM2_z > 1 & !is.na(Erastin_log2))

cat("??? After AIFM2_z > 1 filter:", nrow(filtered), "cell lines remain\n")

# === Step 4: Keep only ModelID, Erastin, and overlapping genes ===
filtered_subset <- filtered %>%
  dplyr::select(ModelID, Erastin_log2, all_of(gene_cols))

# === Step 5: Save filtered dataset ===
write_csv(filtered_subset, "data/processed/ProtExER_Filtered.csv")

cat("??? Filtered dataset saved to ProtExER_Filtered.csv\n")

