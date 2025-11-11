library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)

# ===============================
# Step 1: Load regression results + gene-segment map
# ===============================
reg_results <- read_csv("data/processed/GeneRegressionResults.csv")
gene_segment_map <- read_csv("data/processed/Gene_Segment_Map.csv")

cat("??? Regression results loaded:", nrow(reg_results), "genes\n")
cat("??? Gene-segment map loaded:", nrow(gene_segment_map), "pairs\n")

# ===============================
# Step 2: Merge regression results with segments
# ===============================
reg_results_seg <- reg_results %>%
  left_join(gene_segment_map, by = "Gene")

# ===============================
# Step 3: Summarize per segment
# ===============================
summary_counts <- reg_results_seg %>%
  group_by(SegmentID) %>%
  summarise(
    Positive       = sum(Class == "Positive", na.rm = TRUE),
    Negative       = sum(Class == "Negative", na.rm = TRUE),
    NotSignificant = sum(Class == "NotSignificant", na.rm = TRUE),
    Unknown        = sum(Class == "Unknown" | is.na(Class), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    TotalGenes   = Positive + Negative + NotSignificant + Unknown,
    KnownGenes   = Positive + Negative + NotSignificant,
    Pct_Positive = if_else(KnownGenes > 0, 100 * Positive / KnownGenes, NA_real_),
    Pct_Negative = if_else(KnownGenes > 0, 100 * Negative / KnownGenes, NA_real_)
  ) %>%
  arrange(as.numeric(gsub("seg", "", SegmentID)))   # ??? order by segment number

# Save CSV
write_csv(summary_counts, "data/processed/SegmentSummary.csv")
cat("??? Segment summary saved to SegmentSummary.csv\n")

# ===============================
# Step 4: Barplot (with NotSignificant)
# ===============================
bar_data_all <- summary_counts %>%
  select(SegmentID, Positive, Negative, NotSignificant, KnownGenes) %>%
  pivot_longer(
    cols = Positive:NotSignificant,
    names_to = "Class",
    values_to = "Count"
  ) %>%
  mutate(Fraction = if_else(KnownGenes > 0, Count / KnownGenes, NA_real_))

barplot_all <- ggplot(bar_data_all, aes(x = SegmentID, y = Fraction, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Segment summary (with NotSignificant)",
       y = "Fraction of genes", x = "Segment")

ggsave("figures/PanelA_SegmentSummary_Barplot_ALL.png",
       barplot_all, width = 12, height = 6, dpi = 300)

cat("??? Barplot (with NotSignificant) saved\n")

# ===============================
# Step 5: Barplot (only Significant = Positive/Negative)
# ===============================
bar_data_sig <- summary_counts %>%
  select(SegmentID, Positive, Negative, KnownGenes) %>%
  pivot_longer(
    cols = Positive:Negative,
    names_to = "Class",
    values_to = "Count"
  ) %>%
  mutate(Fraction = if_else(KnownGenes > 0, Count / KnownGenes, NA_real_))

barplot_sig <- ggplot(bar_data_sig, aes(x = SegmentID, y = Fraction, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  labs(title = "Segment summary (only Significant classes)",
       y = "Fraction of genes", x = "Segment")

ggsave("figures/PanelA_SegmentSummary_Barplot_SIG.png",
       barplot_sig, width = 12, height = 6, dpi = 300)

cat("??? Barplot (significant only) saved\n")

# ===============================
# Step 6: Save full summary table (all rows + headers) as PNG
# ===============================
summary_display <- summary_counts %>%
  mutate(
    Pct_Positive = ifelse(!is.na(Pct_Positive), sprintf("%.2f%%", Pct_Positive), NA),
    Pct_Negative = ifelse(!is.na(Pct_Negative), sprintf("%.2f%%", Pct_Negative), NA)
  )

# Create table grob with headers
table_plot <- tableGrob(
  summary_display,
  rows = NULL,
  theme = ttheme_default(
    base_size = 10,
    core = list(fg_params = list(cex = 0.7)),
    colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))
  )
)

# Scale output height to number of rows
nrows <- nrow(summary_display)
png("data/figures/SegmentSummary_Table_FULL.png",
    width = 2000, height = 40 * (nrows + 3), res = 150)
grid.draw(table_plot)
dev.off()

cat("??? Full summary table with headers saved to SegmentSummary_Table_FULL.png\n")
