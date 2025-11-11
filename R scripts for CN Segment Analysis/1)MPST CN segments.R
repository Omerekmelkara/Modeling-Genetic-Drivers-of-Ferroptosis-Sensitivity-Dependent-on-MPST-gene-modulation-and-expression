# Load required packages
library(readr)
library(GenomicRanges)
library(writexl)
library(biomaRt)
library(here)
library(dplyr)   # make sure dplyr is loaded last

# -------------------------------
# 1. Load CN segment data and metadata
# -------------------------------
# Load CN segment data and metadata
cn_segments <- read_csv("data/raw/OmicsCNSegmentsProfileWGS.csv")
omics_profiles <- read_csv("data/raw/OmicsProfiles.csv")
model_metadata <- read_csv("data/raw/model.csv")

# -------------------------------
# 2. Classify CN events (log2-based)
# -------------------------------
cn_segments <- cn_segments %>%
  mutate(event_type = case_when(
    SegmentMean <= -0.5 ~ "Deletion",
    SegmentMean > -0.5 & SegmentMean <= 1 ~ "Gain",
    SegmentMean > 1 ~ "Amplification",
    TRUE ~ "Neutral"
  ))

# -------------------------------
# 3. Map ProfileID ??? ModelID ??? CellLineName
# -------------------------------
cn_segments <- cn_segments %>%
  left_join(dplyr::select(omics_profiles, ProfileID, ModelID), by = "ProfileID") %>%
  left_join(dplyr::select(model_metadata, ModelID, CellLineName, TissueOrigin, OncotreePrimaryDisease),
            by = "ModelID")

# -------------------------------
# 4. Retrieve gene coordinates from Ensembl
# -------------------------------
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
               mart = ensembl)

# -------------------------------
## Remove "chr" prefix if present
cn_segments$Chromosome <- gsub("^chr", "", cn_segments$Chromosome)

# Recreate GenomicRanges
cn_gr <- makeGRangesFromDataFrame(cn_segments,
                                  seqnames.field = "Chromosome",
                                  start.field = "Start",
                                  end.field = "End",
                                  keep.extra.columns = TRUE)

genes_gr <- makeGRangesFromDataFrame(genes,
                                     seqnames.field = "chromosome_name",
                                     start.field = "start_position",
                                     end.field = "end_position",
                                     keep.extra.columns = TRUE)

# Recompute overlaps
overlaps <- findOverlaps(cn_gr, genes_gr)
affected_genes <- data.frame(
  ProfileID = cn_gr$ProfileID[queryHits(overlaps)],
  CellLineName = cn_gr$CellLineName[queryHits(overlaps)],
  TissueOrigin = cn_gr$TissueOrigin[queryHits(overlaps)],
  OncotreePrimaryDisease = cn_gr$OncotreePrimaryDisease[queryHits(overlaps)],
  event_type = cn_gr$event_type[queryHits(overlaps)],
  gene = genes_gr$hgnc_symbol[subjectHits(overlaps)]
)

# Check the first few rows
head(affected_genes)


# -------------------------------
# 7. Summarize CN events by gene and event type
# -------------------------------
summary_table <- affected_genes %>%
  group_by(event_type) %>%
  summarise(
    num_occurrences = n(),
    genes = paste(unique(gene), collapse = ", ")
  )

# Optionally, you can also summarize per cell line
summary_by_cellline <- affected_genes %>%
  group_by(CellLineName, event_type) %>%
  summarise(
    num_occurrences = n(),
    genes = paste(unique(gene), collapse = ", ")
  )

library(dplyr)
library(GenomicRanges)
library(writexl)

# -------------------------------
# 1. Identify CN segments overlapping MPST
# -------------------------------
mpst_gr <- genes_gr[genes_gr$hgnc_symbol == "MPST"]
mpst_overlaps <- findOverlaps(cn_gr, mpst_gr)
mpst_segments <- cn_gr[unique(queryHits(mpst_overlaps))]  # unique segments

# -------------------------------
# 2. Find all genes overlapping these segments
# -------------------------------
all_overlaps <- findOverlaps(mpst_segments, genes_gr)

segment_gene_df <- data.frame(
  segment_index = queryHits(all_overlaps),
  gene = genes_gr$hgnc_symbol[subjectHits(all_overlaps)]
)

# -------------------------------
# 3. Aggregate genes per segment
# -------------------------------
mpst_segment_table <- data.frame(
  ProfileID = mpst_segments$ProfileID,
  CellLineName = mpst_segments$CellLineName,
  TissueOrigin = mpst_segments$TissueOrigin,
  OncotreePrimaryDisease = mpst_segments$OncotreePrimaryDisease,
  Segment_Start = start(mpst_segments),
  Segment_End = end(mpst_segments),
  Segment_Mean = mpst_segments$SegmentMean,
  event_type = mpst_segments$event_type,
  segment_index = seq_along(mpst_segments)  # unique index for joining
)

# Join and collapse gene names per segment
mpst_segment_table <- mpst_segment_table %>%
  left_join(segment_gene_df, by = "segment_index") %>%
  group_by(ProfileID, CellLineName, TissueOrigin, OncotreePrimaryDisease,
           Segment_Start, Segment_End, Segment_Mean, event_type) %>%
  summarise(co_affected_genes = paste(unique(gene), collapse = ", "), .groups = "drop")

# -------------------------------
# 4. Export to CSV
# -------------------------------


# Save processed data
write_csv(mpst_segment_table, "data/processed/MPST_segments_with_genes.csv")


# -------------------------------
# SEGMENT OVERALLS
# -------------------------------

# Compute segment lengths
mpst_segment_table <- mpst_segment_table %>%
  mutate(segment_length = Segment_End - Segment_Start + 1)

# Compute statistics
segment_stats <- mpst_segment_table %>%
  summarise(
    num_segments = n(),
    most_upstream_start = min(Segment_Start),
    most_downstream_end = max(Segment_End),
    shortest_length = min(segment_length),
    longest_length = max(segment_length),
    mean_length = mean(segment_length)
  )

# Print the table
print(segment_stats)
