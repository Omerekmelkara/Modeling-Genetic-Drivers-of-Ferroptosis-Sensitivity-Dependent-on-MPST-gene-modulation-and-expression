library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# === Step 1: Read processed lines ===
raw <- read_lines("data/processed/mpst_cn_filtered.csv")

# Parse the CSV as text
df <- read_csv(
  I(raw),
  col_types = cols(.default = "c"),
  quote = "\""
)

# === Step 2: Extract co_affected_genes from the 9th quoted field ===
# Use grepl to avoid 'sub' returning the entire line on non-matches
pattern <- '^([^,]*,){8}"(.*?)",'
gene_col <- ifelse(
  grepl(pattern, raw[-1]),
  sub(pattern, "\\2", raw[-1]),
  NA_character_
)

df$co_affected_genes <- gene_col

# Quick peek
cat("Preview of co_affected_genes:\n")
print(substr(df$co_affected_genes[1], 1, 200))

# === Step 3: Add segment IDs ===
df <- df %>% mutate(SegmentID = paste0("seg", row_number()))

# === Step 4: Clean + explode into gene???segment pairs ===
gene_segment_map <- df %>%
  mutate(
    co_affected_genes = coalesce(co_affected_genes, ""),
    # remove leading/trailing commas + spaces
    co_affected_genes = str_remove_all(co_affected_genes, "^,\\s*|\\s*,$")
  ) %>%
  separate_rows(co_affected_genes, sep = ",") %>%
  mutate(
    Gene = str_squish(co_affected_genes),
    Gene = str_remove_all(Gene, '^"|"$')   # drop any lingering quotes
  ) %>%
  filter(Gene != "") %>%
  distinct(SegmentID, Gene)

# === Step 5: Save outputs ===
write_csv(gene_segment_map,"data/processed/Gene_Segment_Map.csv")

unique_genes <- distinct(gene_segment_map, Gene)
write_csv(unique_genes, "data/processed/Unique_Genes.csv")

# === Step 6: Report ===
cat("??? Extracted", nrow(gene_segment_map), "gene???segment pairs\n")
cat("??? Unique genes:", nrow(unique_genes), "\n")

# Sanity checks:
gene_segment_map %>% count(SegmentID, sort = TRUE) %>% head(10) %>% print(n = 10)

