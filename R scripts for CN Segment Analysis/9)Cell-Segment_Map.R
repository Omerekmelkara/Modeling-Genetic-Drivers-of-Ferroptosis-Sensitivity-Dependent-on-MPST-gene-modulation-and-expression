library(dplyr)
library(readr)
library(stringr)
library(gridExtra)
library(tidyr)

# === Step 1: Read raw CNV file ===
df <- read_csv(
  "data/processed/mpst_cn_filtered.csv",
  col_types = cols(.default = "c"),
  quote = "\""
)

df <- df %>%
  separate(ProfileID, into = c("ProfileID", "CellLineName", "TissueOrigin", "OncotreePrimaryDisease",
                               "Segment_Start", "Segment_End", "Segment_Mean", "event_type", "co_affected_genes", "ModelID"),
           sep = ",", extra = "merge", fill = "right")

# === Step 2: Add Segment IDs (if not already present) ===
if (!"SegmentID" %in% colnames(df)) {
  df <- df %>% mutate(SegmentID = paste0("seg", row_number()))
}

# === Step 3: Load OmicsProfiles (clean names) ===
omics_info <- read_csv("data/raw/OmicsProfiles.csv") %>%
  mutate(ProfileID = str_trim(ProfileID))

# === Step 4: Drop raw CellLineName from CNV & clean ProfileID ===
df <- df %>%
  select(-CellLineName) %>%
  mutate(ProfileID = str_trim(ProfileID))

# === Step 5: Join to get cleaned CellLineName ===
cell_segment_map <- df %>%
  left_join(
    omics_info %>% transmute(ProfileID, CellLineName = StrippedCellLineName),
    by = "ProfileID"
  ) %>%
  select(CellLineName, SegmentID, event_type)

# === Step 6: Save outputs ===
write_csv(cell_segment_map, "data/processed/Cell_Segment_Map.csv")
cat("??? Saved Cell_Segment_Map.csv with", nrow(cell_segment_map), "rows\n")

# === Step 7: PNG preview with all rows ===
nrows <- nrow(cell_segment_map)
row_height <- 35   # pixels per row (adjust if needed)
png_height <- max(900, nrows * row_height)  # scale height dynamically

png("figures/Cell_Segment_Map.png", width = 1600, height = png_height)
grid.table(cell_segment_map)   # use full dataset, no head()
dev.off()
cat("??? PNG table saved to Cell_Segment_Map.png (", nrows, "rows )\n", sep="")

# === Step 8: Report ===
cat("=== SUMMARY ===\n")
cat("Cell lines:", n_distinct(cell_segment_map$CellLineName), "\n")
cat("Segments:", n_distinct(cell_segment_map$SegmentID), "\n")
cat("Event types:", paste(unique(cell_segment_map$event_type), collapse = ", "), "\n")


