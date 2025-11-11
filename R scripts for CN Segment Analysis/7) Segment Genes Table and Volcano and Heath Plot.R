# =========================
# PANEL B: Volcano Plot of Gene Regression Results + Stats Table
# =========================
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tibble)
library(gridExtra)

# ======================
# 1. Load regression results + segment map
# ======================
reg_results <- read_csv("data/processed/GeneRegressionResults.csv")
gene_segment_map <- read_csv("data/processed/Gene_Segment_Map.csv")

cat("??? Regression results loaded:", nrow(reg_results), "genes\n")
cat("??? Gene-segment map loaded:", nrow(gene_segment_map), "pairs\n")

# ======================
# 2. Filter to genes in segments
# ======================
segment_genes <- unique(gene_segment_map$Gene)

volcano_data <- reg_results %>%
  filter(Gene %in% segment_genes) %>%
  filter(!is.na(Slope), !is.na(Pvalue)) %>%
  mutate(
    PearsonR = sqrt(R2) * sign(Slope)  # reconstruct signed correlation
  )

# ======================
# 3. Counts per class for legend
# ======================
counts <- volcano_data %>%
  count(Class)

class_labels <- setNames(
  paste0(counts$Class, " (", counts$n, ")"),
  counts$Class
)

# ======================
# 4. Top genes for caption
# ======================
top_pos_gene <- volcano_data %>%
  filter(Class == "Positive") %>%
  arrange(Pvalue) %>%
  slice(1) %>%
  pull(Gene)

top_neg_gene <- volcano_data %>%
  filter(Class == "Negative") %>%
  arrange(Pvalue) %>%
  slice(1) %>%
  pull(Gene)

caption_text <- paste0(
  "Top Positive: ", ifelse(length(top_pos_gene)==0, "None", top_pos_gene),
  "   |   Top Negative: ", ifelse(length(top_neg_gene)==0, "None", top_neg_gene)
)

# ======================
# 5. Volcano Plot
# ======================
volcano <- ggplot(volcano_data, aes(x = PearsonR, y = -log10(Pvalue), color = Class)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c("Positive"="red", "Negative"="blue",
               "NotSignificant"="grey70", "Unknown"="black"),
    labels = class_labels
  ) +
  theme_bw() +
  labs(title = "B) Volcano Plot of Segment Genes (Protein vs Erastin, AIFM2 > 1)",
       x = "Signed correlation (Pearson r ~ slope)",
       y = "-log10(p-value)",
       caption = caption_text) +
  theme(legend.position = "top",
        plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))

ggsave("figures/PanelB_SegmentGenes_Volcano.png",
       volcano, width = 9, height = 7, dpi = 300)

cat("??? Volcano plot saved to PanelB_SegmentGenes_Volcano.png\n")

# ======================
# 6. Statistics table for Positive/Negative genes
# ======================
stats_table <- volcano_data %>%
  filter(Class %in% c("Positive", "Negative")) %>%
  arrange(Class, Pvalue) %>%
  select(Gene, Class, Slope, PearsonR, Pvalue, R2, N)

# Save as CSV
write_csv(stats_table, "data/processed/PanelB_PosNeg_Stats.csv")
cat("??? Stats table saved to PanelB_PosNeg_Stats.csv\n")

# Save as PNG (table preview)
stats_plot <- tableGrob(stats_table, rows = NULL, theme = ttheme_default())

ggsave("figures/PanelB_PosNeg_Stats.png",
       stats_plot, width = 12, height = 8, dpi = 300)

cat("??? Stats table with Positive/Negative genes saved to PanelB_PosNeg_Stats.png\n")

# ======================
## ======================
# ======================
# 7. Heatmap of PearsonR (Positive=Red, Negative=Blue, NA=White)
# ======================
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# ====================================
# 7.1. Filter significant genes (Positive/Negative only)
# ====================================
heatmap_data <- reg_results_seg %>%
  filter(Class %in% c("Positive", "Negative")) %>%
  select(Gene, SegmentID, Slope) %>%   # using slope, replace with PearsonR if available
  group_by(Gene, SegmentID) %>%
  summarise(Value = mean(Slope, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SegmentID, values_from = Value, values_fill = NA)

# ====================================
# 7.2. Convert to matrix
# ====================================
heatmap_mat <- as.data.frame(heatmap_data)
rownames(heatmap_mat) <- heatmap_mat$Gene
heatmap_mat <- heatmap_mat[, -1, drop = FALSE]

# Force columns (segments) into numeric order: seg1, seg2, ???
seg_order <- paste0("seg", 1:28)  
seg_order <- intersect(seg_order, colnames(heatmap_mat))  # keep only those present
heatmap_mat <- heatmap_mat[, seg_order, drop = FALSE]

# ====================================
# 7.3. Define color palette (blue ??? white ??? red)
# ====================================
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
value_range <- max(abs(heatmap_mat), na.rm = TRUE)
breaks <- seq(-value_range, value_range, length.out = 101)

# ====================================
# 7.4. Plot & save
# ====================================
png("figures/PanelC_Heatmap.png", width = 1400, height = 1000, res = 150)
pheatmap(
  heatmap_mat,
  color = color_palette,
  breaks = breaks,
  cluster_rows = TRUE,       # cluster genes
  cluster_cols = FALSE,      # keep segments in fixed order
  main = "Heatmap of Significant Genes across Segments",
  fontsize_row = 6, fontsize_col = 8,
  na_col = "white"
)
dev.off()

cat("??? Heatmap saved with segments ordered from seg1 ??? seg28\n")
