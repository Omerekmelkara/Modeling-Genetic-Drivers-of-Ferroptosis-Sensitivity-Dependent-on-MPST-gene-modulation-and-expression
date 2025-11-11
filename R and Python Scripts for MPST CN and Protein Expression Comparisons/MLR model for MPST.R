# ================================
# Full Analysis: MPST, AIFM2, GPX4 & Erastin
# Linear Models + Random Forest + Tables
# ================================

library(tidyverse)
library(broom)
library(ggplot2)
library(cowplot)
library(randomForest)
library(pdp)
library(openxlsx)
library(gridExtra)

# ---------------------------
# 1. Load & preprocess data
# ---------------------------
data <- read.csv("data/processed/combined_protein_erastin_resistance.csv")

colnames(data)[colnames(data) ==
                 "ERASTIN..BRD.BRD.A25004090.001.08.4..log2.fold.change.PRISM.Repurposing.Public.24Q2"] <- "Erastin_log2"

# Z-score normalization + grouping
data <- data %>%
  mutate(
    MPST_z  = as.numeric(scale(MPST)),
    GPX4_z  = as.numeric(scale(GPX4)),
    AIFM2_z = as.numeric(scale(AIFM2)),
    AIFM2_group = case_when(
      AIFM2_z < 0 ~ "Low",
      AIFM2_z >= 0 & AIFM2_z <= 1 ~ "Medium",
      AIFM2_z > 1 ~ "High"
    )
  )

# Theme & colors
custom_theme <- theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        legend.position="right")
aifm2_colors <- c("Low"="#1b9e77","Medium"="#d95f02","High"="#7570b3")

# ---------------------------
# 2. Global + Stratified Linear Models
# ---------------------------

# Global model (All cells)
df_all <- data %>%
  dplyr::select(Erastin_log2, MPST_z, GPX4_z, AIFM2_z) %>%
  drop_na()

mlr_all <- lm(Erastin_log2 ~ MPST_z + GPX4_z + AIFM2_z, data=df_all)

mlr_all_tidy <- broom::tidy(mlr_all, conf.int=TRUE) %>%
  mutate(SampleSize = nrow(df_all),
         AIFM2_group = "All")

# Stratified models
lm_results <- data %>%
  group_by(AIFM2_group) %>%
  group_map(~{
    df <- .x %>%
      dplyr::select(Erastin_log2, MPST_z, GPX4_z, AIFM2_z) %>%
      drop_na()
    this_group <- .y$AIFM2_group
    
    if (nrow(df) >= 4) {
      lm(Erastin_log2 ~ MPST_z + GPX4_z + AIFM2_z, data=df) %>%
        broom::tidy(conf.int=TRUE) %>%
        mutate(SampleSize = nrow(df),
               AIFM2_group = this_group)
    } else {
      tibble() # skip groups too small
    }
  }) %>%
  bind_rows()

# Combine global + stratified
lm_full <- bind_rows(mlr_all_tidy, lm_results)

# ---------------------------
# 3. Scatter + Forest Plots
# ---------------------------

# Scatter plot with adjusted regression lines
pred_data <- data %>%
  drop_na(MPST_z, AIFM2_z, GPX4_z, Erastin_log2, AIFM2_group) %>%
  group_by(AIFM2_group) %>%
  group_modify(~{
    mod <- lm(Erastin_log2 ~ MPST_z + GPX4_z + AIFM2_z, data=.x)
    .x %>% mutate(Predicted = predict(mod, newdata=.x))
  }) %>%
  ungroup()

scatter_mlr_plot <- ggplot(pred_data, aes(x=MPST_z, y=Erastin_log2, color=AIFM2_group)) +
  geom_point(alpha=0.6, size=2) +
  geom_line(aes(y=Predicted), size=1.2) +
  facet_wrap(~AIFM2_group) +
  scale_color_manual(values=aifm2_colors) +
  labs(x="MPST (z-score)", y="Erastin sensitivity (log2)", color="AIFM2 group") +
  custom_theme +
  ggtitle("Stratified MLR: Adjusted MPST effect on Erastin")

# Forest plot ??? MPST only
forest_mpst <- lm_full %>%
  filter(term=="MPST_z") %>%
  ggplot(aes(x=estimate, y=AIFM2_group, color=AIFM2_group)) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0.25, size=1) +
  geom_vline(xintercept=0, linetype="dashed", color="gray40") +
  scale_color_manual(values=c("All"="black", aifm2_colors)) +
  labs(x="Slope (MPST effect)", y="Group", color="Group") +
  custom_theme +
  ggtitle("Forest plot of MPST effect (Global + Stratified)")

# Combined figure
fig1 <- plot_grid(scatter_mlr_plot, forest_mpst, labels=c("A","B"), ncol=2, rel_widths=c(1.4,1))

# Forest plot ??? all predictors
forest_all <- lm_full %>%
  filter(term!="(Intercept)") %>%
  ggplot(aes(x=estimate, y=term, color=AIFM2_group)) +
  geom_point(size=3, position=position_dodge(width=0.6)) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high),
                 position=position_dodge(width=0.6), height=0.2) +
  geom_vline(xintercept=0, linetype="dashed", color="gray40") +
  scale_color_manual(values=c("All"="black", aifm2_colors)) +
  labs(x="Effect size (slope)", y="Predictor", color="Group") +
  custom_theme +
  ggtitle("Forest plots: predictor effects (Global + Stratified)")

# ---------------------------
# 4. Random Forest Importance
# ---------------------------
make_rf_importance <- function(df, title, fillcol) {
  rf_model <- randomForest(Erastin_log2 ~ ., data=df, importance=TRUE, ntree=1000)
  imp_df <- importance(rf_model, type=1) %>%
    as.data.frame() %>%
    rownames_to_column("Gene")
  colnames(imp_df)[2] <- "Importance"
  imp_df <- imp_df %>% arrange(desc(Importance))
  list(
    model = rf_model,
    data = imp_df,
    plot = ggplot(imp_df, aes(x=reorder(Gene, Importance), y=Importance)) +
      geom_col(fill=fillcol) +
      coord_flip() +
      labs(x="Gene", y="Importance (MDA)") +
      custom_theme + ggtitle(title)
  )
}

df_low   <- data %>% filter(AIFM2_group=="Low") %>% dplyr::select(Erastin_log2, MPST_z, AIFM2_z, GPX4_z) %>% drop_na()
df_med   <- data %>% filter(AIFM2_group=="Medium") %>% dplyr::select(Erastin_log2, MPST_z, AIFM2_z, GPX4_z) %>% drop_na()
df_high  <- data %>% filter(AIFM2_group=="High") %>% dplyr::select(Erastin_log2, MPST_z, AIFM2_z, GPX4_z) %>% drop_na()

rf_all   <- make_rf_importance(df_all,   "RF importance (All cells)",    "steelblue")
rf_low   <- make_rf_importance(df_low,   "RF importance (Low AIFM2)",    "#1b9e77")
rf_med   <- make_rf_importance(df_med,   "RF importance (Medium AIFM2)", "#d95f02")
rf_high  <- make_rf_importance(df_high,  "RF importance (High AIFM2)",   "#7570b3")

# ---------------------------
# 5. PDP Panels (All cells model)
# ---------------------------
df_all_plot <- data %>%
  dplyr::select(Erastin_log2, MPST_z, AIFM2_z, GPX4_z, AIFM2_group) %>%
  drop_na()

pdp_mpst <- partial(rf_all$model, pred.var="MPST_z", train=df_all, grid.resolution=50) %>%
  as_tibble() %>%
  ggplot(aes(x=MPST_z, y=yhat)) +
  geom_point(data=df_all_plot, aes(x=MPST_z, y=Erastin_log2, color=AIFM2_group),
             inherit.aes=FALSE, alpha=0.6, size=2) +
  geom_line(color="steelblue", size=1.2) +
  scale_color_manual(values=aifm2_colors) +
  labs(x="MPST (z-score)", y="Erastin sensitivity (log2)",
       title="Panel A: RF PDP for MPST_z", color="AIFM2 group") +
  custom_theme

pdp_aifm2 <- partial(rf_all$model, pred.var="AIFM2_z", train=df_all, grid.resolution=50) %>%
  as_tibble() %>%
  ggplot(aes(x=AIFM2_z, y=yhat)) +
  geom_point(data=df_all_plot, aes(x=AIFM2_z, y=Erastin_log2, color=AIFM2_group),
             inherit.aes=FALSE, alpha=0.6, size=2) +
  geom_line(color="darkorange", size=1.2) +
  scale_color_manual(values=aifm2_colors) +
  labs(x="AIFM2 (z-score)", y="Erastin sensitivity (log2)",
       title="Panel B: RF PDP for AIFM2_z", color="AIFM2 group") +
  custom_theme

pdp_gpx4 <- partial(rf_all$model, pred.var="GPX4_z", train=df_all, grid.resolution=50) %>%
  as_tibble() %>%
  ggplot(aes(x=GPX4_z, y=yhat)) +
  geom_point(data=df_all_plot, aes(x=GPX4_z, y=Erastin_log2, color=AIFM2_group),
             inherit.aes=FALSE, alpha=0.6, size=2) +
  geom_line(color="purple", size=1.2) +
  scale_color_manual(values=aifm2_colors) +
  labs(x="GPX4 (z-score)", y="Erastin sensitivity (log2)",
       title="Panel C: RF PDP for GPX4_z", color="AIFM2 group") +
  custom_theme

interaction_pdp <- partial(rf_all$model, pred.var=c("MPST_z","AIFM2_z"), train=df_all,
                           grid.resolution=30, chull=TRUE)
interaction_df <- as_tibble(interaction_pdp)

interaction_plot <- ggplot(interaction_df, aes(x=MPST_z, y=AIFM2_z, fill=yhat)) +
  geom_tile() +
  geom_contour(aes(z=yhat), color="black", alpha=0.5) +
  scale_fill_viridis_c(option="plasma") +
  labs(x="MPST (z-score)", y="AIFM2 (z-score)", fill="Pred. Erastin",
       title="Panel D: RF 2D PDP (MPST ?? AIFM2)") +
  custom_theme

# ---------------------------
# 6. Save outputs
# ---------------------------
output_dir <- "C:/Users/Ekmel/Desktop/Plots"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# Linear models
ggsave(file.path(output_dir,"Fig1_Stratifed_LM.png"), fig1, width=13, height=5, dpi=300)
ggsave(file.path(output_dir,"Fig1b_Forest_AllPredictors.png"), forest_all, width=9, height=5.5, dpi=300)

# Tables for LM
lm_table_mpst <- lm_full %>%
  filter(term=="MPST_z") %>%
  select(AIFM2_group, SampleSize, estimate, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x,3)))
ggsave(file.path(output_dir,"Fig1_MPST_Effect_Table.png"),
       tableGrob(lm_table_mpst, rows=NULL), width=8, height=3, dpi=300)

lm_table_all <- lm_full %>%
  filter(term!="(Intercept)") %>%
  select(AIFM2_group, term, SampleSize, estimate, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x,3)))
ggsave(file.path(output_dir,"Fig1b_AllPredictors_Table.png"),
       tableGrob(lm_table_all, rows=NULL), width=10, height=5, dpi=300)

write.xlsx(lm_full, file.path(output_dir,"Stratified_LM_GlobalPlusSubgroups.xlsx"))

# Random Forest plots + tables
plots_rf <- list(all=rf_all, low=rf_low, med=rf_med, high=rf_high)
for (nm in names(plots_rf)) {
  ggsave(file.path(output_dir,paste0("Fig2_RF_Importance_",nm,".png")),
         plots_rf[[nm]]$plot, width=7, height=5, dpi=300)
  ggsave(file.path(output_dir,paste0("Fig2_RF_Importance_",nm,"_Table.png")),
         tableGrob(plots_rf[[nm]]$data, rows=NULL), width=6, height=3, dpi=300)
}

# PDP Panels
ggsave(file.path(output_dir,"PanelA_RF_PDP_MPST.png"), pdp_mpst, width=6, height=5, dpi=300)
ggsave(file.path(output_dir,"PanelB_RF_PDP_AIFM2.png"), pdp_aifm2, width=6, height=5, dpi=300)
ggsave(file.path(output_dir,"PanelC_RF_PDP_GPX4.png"), pdp_gpx4, width=6, height=5, dpi=300)
ggsave(file.path(output_dir,"PanelD_RF_PDP_Interaction.png"), interaction_plot, width=7, height=6, dpi=300)

message("??? All figures and tables saved to: ", output_dir)


# ---------------------------
# Ensure "All" group exists
# ---------------------------
data_with_all <- data %>%
  mutate(AIFM2_group = ifelse(is.na(AIFM2_group), "All", AIFM2_group))

# Add global "All" group explicitly
data_with_all <- bind_rows(
  data_with_all,
  data %>% mutate(AIFM2_group = "All")
)

# ---------------------------
# Function for stratified MLR scatter plots
# ---------------------------
make_stratified_mlr_plot <- function(df, predictor, title, color_line="steelblue") {
  ggplot(df, aes_string(x=predictor, y="Erastin_log2", color="AIFM2_group")) +
    geom_point(alpha=0.6, size=2) +
    geom_smooth(method="lm", se=TRUE, size=1.2, color=color_line) +
    facet_wrap(~AIFM2_group) +
    scale_color_manual(values=c(aifm2_colors, All="black")) +
    labs(x=paste0(predictor, " (z-score)"), 
         y="Erastin sensitivity (log2)", 
         color="AIFM2 group") +
    custom_theme +
    ggtitle(title)
}

# ---------------------------
# Generate plots for each predictor
# ---------------------------
mlr_plot_mpst  <- make_stratified_mlr_plot(data_with_all, "MPST_z",  "Stratified MLR: MPST vs Erastin",  "steelblue")
mlr_plot_aifm2 <- make_stratified_mlr_plot(data_with_all, "AIFM2_z", "Stratified MLR: AIFM2 vs Erastin", "darkorange")
mlr_plot_gpx4  <- make_stratified_mlr_plot(data_with_all, "GPX4_z",  "Stratified MLR: GPX4 vs Erastin",  "purple")

# ---------------------------
# Save figures
# ---------------------------
ggsave(file.path(output_dir,"StratifiedMLR_MPST.png"),  mlr_plot_mpst,  width=7, height=5, dpi=300)
ggsave(file.path(output_dir,"StratifiedMLR_AIFM2.png"), mlr_plot_aifm2, width=7, height=5, dpi=300)
ggsave(file.path(output_dir,"StratifiedMLR_GPX4.png"),  mlr_plot_gpx4,  width=7, height=5, dpi=300)
