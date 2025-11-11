library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(gridExtra)
library(grid)

# ===============================
# Step 1: Load data
# ===============================
prot_er <- read_csv("data/processed/ProtExER_Filtered.csv")
gene_segment_map <- read_csv("data/processed/Gene_Segment_Map.csv")
cn_data <- read_csv("data/processed/mpst_cn_filtered.csv") 

# ===============================
# Step 2: Long format proteomics
# ===============================
prot_er_long <- prot_er %>%
  select(ModelID, ends_with("_z")) %>%
  pivot_longer(
    cols = ends_with("_z"),
    names_to = "Gene",
    values_to = "ProteinZ"
  ) %>%
  mutate(Gene = str_remove(Gene, "_z$"))

# ===============================
# Step 3: Map genes ??? segments
# ===============================
prot_seg <- prot_er_long %>%
  inner_join(gene_segment_map, by = "Gene")

# ===============================
# Step 4: Add CellLineName lookup
# ===============================
cell_lookup <- cn_data %>%
  select(ModelID, CellLineName) %>%
  distinct()

prot_seg <- prot_seg %>%
  left_join(cell_lookup, by = "ModelID")

# ===============================
# Step 5: Collapse to Segment ??? (Model/CellLine)
# ===============================
segment_cells <- prot_seg %>%
  distinct(SegmentID, ModelID, CellLineName) %>%
  group_by(SegmentID) %>%
  summarise(
    Models    = paste(unique(ModelID), collapse = "; "),
    CellLines = paste(unique(CellLineName), collapse = "; "),
    Nmodels   = n_distinct(ModelID),
    .groups = "drop"
  )

# ===============================
# Step 6: Save
# ===============================
write_csv(segment_cells, "data/processed/Segment_to_CellPairs.csv")
cat("??? Segment-to-cell mapping saved with models + cell lines\n")

# ===============================
# Step 7: Big PNG table (Cells only, no Model IDs)
# ===============================
seg_cells_tbl <- segment_cells %>%
  mutate(
    Segment   = SegmentID,
    `# Models` = Nmodels,
    CellLines = str_squish(str_wrap(CellLines, width = 120))
  ) %>%
  select(Segment, `# Models`, CellLines) %>%
  arrange(as.numeric(gsub("\\D", "", Segment)))  # seg1???segN order

tab <- tableGrob(
  seg_cells_tbl,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(fg_params = list(fontsize = 8, lineheight = 1.1),
                padding = unit(c(6, 6), "pt")),
    colhead = list(fg_params = list(fontsize = 10, fontface = "bold"),
                   bg_params = list(fill = "#F2F2F2"))
  )
)

n_rows <- nrow(seg_cells_tbl)
png_width_in  <- 20
png_height_in <- max(8, min(60, 0.35 * n_rows + 3))

png("figures/Segment_to_CellPairs_CellsOnly.png",
    width = png_width_in, height = png_height_in, units = "in", res = 200)
grid.newpage(); grid.draw(tab)
dev.off()

cat("??? PNG saved: Segment_to_CellPairs_CellsOnly.png\n")
