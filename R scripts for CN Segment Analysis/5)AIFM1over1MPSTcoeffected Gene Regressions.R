library(dplyr)
library(broom)
library(purrr)
library(readr)
library(stringr)

# === Step 1: Load the filtered data ===
prot_er <- read_csv("data/processed/depmap data/ProtExER_Filtered.csv")

# Get list of protein z-score columns (they should end with "_z")
protein_cols <- grep("_z$", colnames(prot_er), value = TRUE)

# Set threshold
min_n <- 30

# === Step 2: Run regression per gene ===
reg_results <- map_dfr(protein_cols, function(gene_col) {
  # Build safe formula using backticks
  f <- as.formula(paste("Erastin_log2 ~ `", gene_col, "`", sep = ""))
  
  # Prepare data
  dat <- prot_er %>%
    select(all_of(c("Erastin_log2", gene_col))) %>%
    filter(!is.na(Erastin_log2) & !is.na(.data[[gene_col]]))
  
  n_obs <- nrow(dat)
  
  if (n_obs >= min_n) {
    fit <- tryCatch(lm(f, data = dat), error = function(e) NULL)
    
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit) %>% filter(term == paste0("`", gene_col, "`") | term == gene_col)
      glance_fit <- glance(fit)
      
      tibble(
        Gene = sub("_z$", "", gene_col),
        Slope = tidy_fit$estimate,
        Pvalue = tidy_fit$p.value,
        R2 = glance_fit$r.squared,
        N = n_obs
      )
    } else {
      tibble(Gene = sub("_z$", "", gene_col), 
             Slope = NA, Pvalue = NA, R2 = NA, N = n_obs)
    }
  } else {
    tibble(Gene = sub("_z$", "", gene_col), 
           Slope = NA, Pvalue = NA, R2 = NA, N = n_obs)
  }
})


# === Step 3: Classification ===
reg_results <- reg_results %>%
  mutate(
    Class = case_when(
      is.na(Pvalue) | N < min_n ~ "Unknown",
      Pvalue < 0.05 & Slope > 0 ~ "Positive",
      Pvalue < 0.05 & Slope < 0 ~ "Negative",
      TRUE ~ "NotSignificant"
    )
  )

# === Step 4: Save results ===
write_csv(reg_results,"data/processed/GeneRegressionResults.csv")

cat("??? Finished regression on", length(protein_cols), "genes\n")
cat("??? Passed threshold:", sum(reg_results$N >= min_n), "genes\n")
cat("??? Results saved to GeneRegressionResults.csv\n")
