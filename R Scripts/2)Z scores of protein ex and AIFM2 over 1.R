library(here)
library(readr)
library(tibble)
library(dplyr)   # load last

# ==============================
# 1. Load Gygi proteomics
# ==============================
raw <- read_csv("data/raw/gygi_with_symbols.csv", col_names = FALSE)

# First row = protein names (symbols)
protein_names <- raw[1, ] %>% unlist() %>% as.character()

# Assign as colnames
colnames(raw) <- protein_names

# Drop the first 3 annotation rows, keep from row 4 = data
df <- raw[-c(1:3), ]

# Ensure first column is ModelID
colnames(df)[1] <- "ModelID"

# Clean duplicate/empty protein names
clean_names <- colnames(df)
clean_names[is.na(clean_names) | clean_names == ""] <- paste0("Unknown_", seq_len(sum(is.na(clean_names) | clean_names == "")))
clean_names <- make.unique(clean_names)   # add .1, .2 etc. if duplicated
colnames(df) <- clean_names

# Convert numeric columns properly
df[ , -1] <- lapply(df[ , -1], as.numeric)

cat("??? Gygi dimensions:", dim(df), "\n")
print(head(df[, 1:5]))


# ==============================
# 2. Compute z-scores
# ==============================
df_raw <- df

z_scores <- as.data.frame(scale(df_raw[ , -1], center = TRUE, scale = TRUE))
colnames(z_scores) <- paste0(colnames(df_raw)[-1], "_z")

# add ModelID back to z-scores
z_scores <- cbind(ModelID = df_raw$ModelID, z_scores)

# --- merge raw + z ---
df_all <- df_raw %>%
  left_join(z_scores, by = "ModelID")

cat("??? Combined proteomics dimensions:", dim(df_all), "\n")


# ==============================
# 3. Load Erastin data
# ==============================
erastin <- read_csv("data/raw/ERASTIN (BRD_BRD-A25004090-001-08-4) log2 fold change.csv")

# rename erastin column
colnames(erastin)[colnames(erastin) ==
                    "ERASTIN (BRD:BRD-A25004090-001-08-4) log2 fold change PRISM Repurposing Public 24Q2"] <- "Erastin_log2"

cat("??? Erastin dimensions:", dim(erastin), "\n")


# ==============================
# 4. Merge proteomics + Erastin
# ==============================
merged <- df_all %>%
  inner_join(erastin %>% dplyr::select(model, Erastin_log2),
             by = c("ModelID" = "model"))

cat("??? Merged dataset:", nrow(merged), "samples with both proteomics + Erastin\n")


# ==============================
# 5. Subset high-AIFM2 (z > 1)
# ==============================
if(!"AIFM2_z" %in% colnames(merged)){
  stop("??? Column AIFM2_z not found. Check if Gygi file has 'AIFM2' or alias 'FSP1'.")
}

expr_high <- merged %>%
  filter(AIFM2_z > 1 & !is.na(Erastin_log2))

cat("??? High AIFM2 samples with Erastin data:", nrow(expr_high), "\n")

# preview key columns
print(head(expr_high[, c("ModelID", "AIFM2", "AIFM2_z", "Erastin_log2")]))

# Save merged proteomics + Erastin dataset
write_csv(merged, "data/processed/Merged_Proteomics_Erastin.csv")

cat("??? Merged dataset saved to Merged_Proteomics_Erastin.csv\n")


# ==============================
# 5. Load MPST CN data
# ==============================
mpst_cn <- read_csv("data/processed/MPST_segments_with_genes.csv")

cat("??? MPST CN dimensions:", dim(mpst_cn), "\n")
head(mpst_cn)

# ==============================
# 6. Load OmicsProfiles (ProfileID ??? ModelID mapping)
# ==============================
omics_profiles <- read_csv("data/raw/OmicsProfiles.csv")

cat("??? OmicsProfiles dimensions:", dim(omics_profiles), "\n")
head(omics_profiles)

# ==============================
# 7. Map MPST CN to ModelID
# ==============================
mpst_cn_mapped <- mpst_cn %>%
  inner_join(omics_profiles %>% dplyr::select(ProfileID, ModelID),
             by = "ProfileID")

cat("??? Mapped CN entries with ModelID:", nrow(mpst_cn_mapped), "\n")

# ==============================
# 8. Filter CN data to high-AIFM2 cells
# ==============================
mpst_cn_filtered <- mpst_cn_mapped %>%
  filter(ModelID %in% expr_high$ModelID)

cat("??? Final MPST CN rows for high-AIFM2 cells:", nrow(mpst_cn_filtered), "\n")

# preview a few rows
head(mpst_cn_filtered)

# ==============================
# 9. Save filtered CN events
# ==============================
write_csv(mpst_cn_filtered,"data/processed/mpst_cn_filtered.csv")

cat("??? Saved mpst_cn_filtered.csv with", nrow(mpst_cn_filtered), "rows\n")


