library(dplyr)
library(broom)
library(knitr)
library(kableExtra)
library(ggplot2)

# ---------------------------
# Parameters
# ---------------------------
window_size <- 40
step_size <- 5

# ---------------------------
# Load dataset
# ---------------------------
data <- read.csv("data/processed/Merged_CN_Protein_complete.csv")

# Ensure z-scores exist
if(!"AIFM2_z" %in% colnames(data)) {
  data$AIFM2_z <- scale(data$AIFM2_protex, center=TRUE, scale=TRUE)[,1]
}
if(!"MPST_z" %in% colnames(data)) {
  data$MPST_z <- scale(data$MPST_protex, center=TRUE, scale=TRUE)[,1]
}

# Sort by AIFM2_z
data_sorted <- data %>% arrange(AIFM2_z)

# ---------------------------
# Function for sliding window regression
# ---------------------------
window_regression <- function(df, predictor, outcome, x_align, window_size, step_size) {
  results <- data.frame()
  n <- nrow(df)
  
  for(start_idx in seq(1, n - window_size + 1, by = step_size)) {
    window_data <- df[start_idx:(start_idx + window_size - 1), ]
    mod <- lm(as.formula(paste(outcome, "~", predictor)), data = window_data)
    tidy_mod <- tidy(mod)
    
    slope <- tidy_mod$estimate[tidy_mod$term == predictor]
    se    <- tidy_mod$std.error[tidy_mod$term == predictor]
    pval  <- tidy_mod$p.value[tidy_mod$term == predictor]
    
    mean_align <- mean(window_data[[x_align]], na.rm=TRUE)
    n_samples <- nrow(window_data)
    
    results <- rbind(results, data.frame(
      AIFM2_z = mean_align,
      slope = slope,
      se = se,
      p_value = pval,
      n_samples = n_samples,
      predictor = predictor
    ))
  }
  
  return(results)
}

# ---------------------------
# Run sliding regressions
# ---------------------------
protein_window <- window_regression(data_sorted, "MPST_z", "Erastin_log2", "AIFM2_z", window_size, step_size)
cn_window      <- window_regression(data_sorted, "MPST_CN", "Erastin_log2", "AIFM2_z", window_size, step_size)

# ---------------------------
# Plot Protein vs Erastin
# ---------------------------
plot_protein <- ggplot(protein_window, aes(x=AIFM2_z, y=slope)) +
  geom_line(color="steelblue", size=1.2) +
  geom_ribbon(aes(ymin=slope-se, ymax=slope+se), alpha=0.2, fill="steelblue") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  # highlight significant points
  geom_point(data = subset(protein_window, p_value < 0.05),
             aes(x=AIFM2_z, y=slope),
             color="red", size=2.5) +
  labs(x="Mean AIFM2 (z-score, windowed)", 
       y="Slope of MPST_z vs Erastin", 
       title="Sliding-Window Regression: MPST Protein Effect",
       subtitle="Red dots = significant slopes (p < 0.05)") +
  theme_minimal()

# ---------------------------
# Plot CN vs Erastin
# ---------------------------
plot_cn <- ggplot(cn_window, aes(x=AIFM2_z, y=slope)) +
  geom_line(color="darkorange", size=1.2) +
  geom_ribbon(aes(ymin=slope-se, ymax=slope+se), alpha=0.2, fill="darkorange") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  # highlight significant points
  geom_point(data = subset(cn_window, p_value < 0.05),
             aes(x=AIFM2_z, y=slope),
             color="red", size=2.5) +
  labs(x="Mean AIFM2 (z-score, windowed)", 
       y="Slope of MPST_CN vs Erastin", 
       title="Sliding-Window Regression: MPST Copy Number Effect",
       subtitle="Red dots = significant slopes (p < 0.05)") +
  theme_minimal()

# Save plots
ggsave("Protein_vs_Erastin_Windowed.png", plot_protein, width=7, height=5, dpi=300)
ggsave("CN_vs_Erastin_Windowed.png", plot_cn, width=7, height=5, dpi=300)
