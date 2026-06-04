rm(list=ls())

################################################################################
## Master script: Recompute Simulation Results & Generate Figures
##
## This script runs the simulation studies and regenerates the corresponding
## figures/tables block by block, strictly following the flow of the main
## manuscript first, followed by the Supplementary Material.
##
## [REVIEWER INSTRUCTIONS]
## - Full recomputation with default settings (M_pred=500, M_test=10000) may
##   take several days.
## - For a fast reproducibility check, please decrease the 'M_pred' and 'M_test'
##   parameters in the "USER SETTINGS" section below (e.g., M_pred=10, M_test=100).
## - You can run each section individually to instantly view its corresponding figure.
################################################################################

# ==============================================================================
# 0. Setup & Utility
# ==============================================================================
get_current_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }
  return(getwd())
}

current_dir <- get_current_dir()
cat("Script directory:", current_dir, "\n")
setwd(current_dir)

# Load simulation and plotting functions
source("./simulation/simulations.R")
source("./simulation/generate_figs_tabs.R")


# ==============================================================================
# 1. USER SETTINGS (Adjustable for fast reproduction)
# ==============================================================================
# Original manuscript values: M_pred = 500, M_test = 10000.
# Change these values to run a smaller subset for quick verification.
M_pred <- 500   # Number of simulations for prediction & tuning studies
M_test <- 10000 # Number of simulations for unbiasedness tests


# ==============================================================================
# PART A: Main Manuscript Simulations & Figures
# ==============================================================================

# ------------------------------------------------------------------------------
# Section 3.1.2: Recovering the correct tree structure (Manuscript Figure 2)
# ------------------------------------------------------------------------------
cat("\n--- Running Section 3.1.2: Recovering the correct tree structure ---\n")
sim_set_recovery <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                         cens_rates = c(20, 50), ksi = 500, model = "tree")
for(Dist in c("WI", "WD", "Lgn", "Bat")) {
  sim_set_recovery$Dist <- Dist
  simulate_LBRC_tree_methods("ANOVA study", sim_set_recovery, current_dir)
}
# Output: Manuscript Figure 2
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_2_compare_recovery_rate", results_dir)


# ------------------------------------------------------------------------------
# Section 3.1.3: Prediction accuracy against LTRC-CIT (Manuscript Figure 3)
# ------------------------------------------------------------------------------
cat("\n--- Running Section 3.1.3: Prediction accuracy against LTRC-CIT ---\n")
sim_set_anova <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                      cens_rates = c(20, 50), ksi = 500)
for(Dist in c("WI", "WD", "Lgn", "Bat")) {
  sim_set_anova$Dist <- Dist
  if(Dist %in% c("WI", "WD")){
    for(model in c("tree", "linear", "nonlinear", "interaction")) {
      sim_set_anova$model <- model
      simulate_LBRC_tree_methods("ANOVA study", sim_set_anova, current_dir)
    }
  }else{ # Dist %in% c("Lgn", "Bat")
    sim_set_anova$model <- "tree"
    simulate_LBRC_tree_methods("ANOVA study", sim_set_anova, current_dir)
  }
}
# Output: Manuscript Figure 3
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_3_LBRCtrees_vs_LTRCtree_ANOVA", results_dir)


# ------------------------------------------------------------------------------
# Section 3.2.1: Regulating the construction of trees in forests
#                - WI Setting (Manuscript Figure 4)
# ------------------------------------------------------------------------------
cat("\n--- Running Section 3.2.1: Regulating the construction of trees in forests (WI Setting) ---\n")
sim_set_tune_main <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                          cens_rates = 20, ksi = 500, tune.metric = "brier")
for(model in c("tree", "linear", "nonlinear", "interaction")) {
  sim_set_tune_main$model <- model
  sim_set_tune_main$Dist <- "WI"
  simulate_LBRC_tree_methods("OOB tuning validation", sim_set_tune_main, current_dir)
}
# Output: Manuscript Figure 4
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_4_WI_LBRCforests_OOB_tuning_brier", results_dir)


# ------------------------------------------------------------------------------
# Section 3.2.2: Prediction accuracy across methods
#                - WI Setting, 20% Censoring (Manuscript Figure 5)
# ------------------------------------------------------------------------------
cat("\n--- Running Section 3.2.2: Prediction accuracy across methods (WI Setting, 20% Censoring) ---\n")
sim_set_pred_main <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                          cens_rates = 20, ksi = 500)
for(model in c("tree", "linear", "nonlinear", "interaction")) {
  sim_set_pred_main$model <- model
  sim_set_pred_main$Dist <- "WI"
  simulate_LBRC_tree_methods("model prediction", sim_set_pred_main, current_dir)
}
# Output: Manuscript Figure 5
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_5_WI_20_compare_prediction_accuracy", results_dir)



# ==============================================================================
# PART B: Supplementary Material Simulations & Figures
# ==============================================================================

# ------------------------------------------------------------------------------
# Section B: Test of unbiasedness of variable selection (Supp Figure S1)
# ------------------------------------------------------------------------------
cat("\n--- Running Section B: Test of unbiasedness of variable selection ---\n")
sim_set_unbias <- list(M = M_test, ns = 200L, cens_rates = c(20, 50), ksi = 500)
for(Dist in c("WI", "WD", "Lgn")) {
  sim_set_unbias$Dist <- Dist
  simulate_LBRC_tree_methods("test unbiasedness", sim_set_unbias, current_dir)
}
# Output: Supplementary Figure S1
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_S1_test_unbiasedness_LBRCtrees", results_dir)


# ------------------------------------------------------------------------------
# Section C: Additional results on regulating the construction of trees in forests
#            (Supp Figures S2-S6)
# ------------------------------------------------------------------------------
cat("\n--- Running Section C.2.1: Results based on IBS ---\n")
sim_set_tune_supp <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                          cens_rates = 20, ksi = 500, tune.metric = "brier")
for(model in c("tree", "linear", "nonlinear", "interaction")) {
  sim_set_tune_supp$model <- model
  Dists <- if(model == "tree") c("WD", "Lgn", "Bat") else "WD"
  for(Dist in Dists) {
    sim_set_tune_supp$Dist <- Dist
    simulate_LBRC_tree_methods("OOB tuning validation", sim_set_tune_supp, current_dir)
  }
}
# Output: Supplementary Figures S2-S3
results_dir <- file.path(current_dir, "results_intermediate")
for (fig in c("figure_S2_WD_LBRCforests_OOB_tuning_brier",
              "figure_S3_LgnBat_LBRCforests_OOB_tuning_brier")) {
  generate_figures_tables(fig, results_dir)
}

cat("\n--- Running Section C.2.2: Results based on C-index ---\n")
sim_set_tune_supp <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                          cens_rates = 20, ksi = 500, tune.metric = "cindex")
for(model in c("tree", "linear", "nonlinear", "interaction")) {
  sim_set_tune_supp$model <- model
  Dists <- if(model == "tree") c("WI", "WD", "Lgn", "Bat") else c("WI", "WD")
  for(Dist in Dists) {
    sim_set_tune_supp$Dist <- Dist
    simulate_LBRC_tree_methods("OOB tuning validation", sim_set_tune_supp, current_dir)
  }
}
# Output: Supplementary Figures S4-S6
results_dir <- file.path(current_dir, "results_intermediate")
for (fig in c("figure_S4_WI_LBRCforests_OOB_tuning_cindex",
              "figure_S5_WD_LBRCforests_OOB_tuning_cindex",
              "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")) {
  generate_figures_tables(fig, results_dir)
}


# ------------------------------------------------------------------------------
# Section D: Additional results on prediction accuracy across methods
#            (Supp Figures S7-S10)
# ------------------------------------------------------------------------------
cat("\n--- Running Section D: Additional prediction accuracy across methods ---\n")
sim_set_pred_supp <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400), ksi = 500)
for(dist in c("WI", "WD", "Lgn", "Bat")){
  if(dist == "WI"){
    for(model in c("tree", "linear", "nonlinear", "interaction")) {
      sim_set_pred_supp$model <- model
      sim_set_pred_supp$Dist <- dist
      sim_set_pred_supp$cens_rates <- 50
      simulate_LBRC_tree_methods("model prediction", sim_set_pred_supp, current_dir)
    }
  }else if(dist == "WD"){
    for(model in c("tree", "linear", "nonlinear", "interaction")) {
      sim_set_pred_supp$model <- model
      sim_set_pred_supp$Dist <- dist
      sim_set_pred_supp$cens_rates <- c(20, 50)
      simulate_LBRC_tree_methods("model prediction", sim_set_pred_supp, current_dir)
    }
  }else{ # dist %in% c("Lgn", "Bat")
    sim_set_pred_supp$model <- "tree"
    sim_set_pred_supp$Dist <- dist
    sim_set_pred_supp$cens_rates <- c(20, 50)
    simulate_LBRC_tree_methods("model prediction", sim_set_pred_supp, current_dir)
  }
}
# Output: Supplementary Figures S7-S10
results_dir <- file.path(current_dir, "results_intermediate")
for (fig in c("figure_S7_WI_50_compare_prediction_accuracy",
              "figure_S8_WD_20_compare_prediction_accuracy",
              "figure_S9_WD_50_compare_prediction_accuracy",
              "figure_S10_LgnBat_2050_compare_prediction_accuracy")) {
  generate_figures_tables(fig, results_dir)
}


# ------------------------------------------------------------------------------
# Section E: Sensitivity Analysis (Supp Figures S12-S13)
# ------------------------------------------------------------------------------
sim_set_sens <- list(ns = 200L, cens_rates = 20, Dist = "WI")
cat("\n--- Running Section E.1: Unbiasedness of variable selection ---\n")
for(scenario in c("unbias_texpt", "unbias_covd")) {
  sim_set_sens$scenario <- scenario
  sim_set_sens$M <- M_test
  sim_set_sens$tau <- qweibull(0.9999, 2, 3) + 1
  for(mu in c(0.1, 0.2, 0.5)) {
    sim_set_sens$rho <- sim_set_sens$tau * mu
    simulate_LBRC_tree_methods("sensitivity analysis", sim_set_sens, current_dir)
  }
}
# Output: Supplementary Figure S12
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_S11_sensitivity_analysis_unbiasedness", results_dir)

cat("\n--- Running Section E.2: Tree recovery and prediction accuracy ---\n")
for(scenario in c("tree_texpt", "tree_covd", "nlin_texpt", "nlin_covd")) {
  sim_set_sens$scenario <- scenario
  sim_set_sens$M <- M_pred
  sim_set_sens$cov_set_num <- 10
  if(scenario %in% c("tree_texpt", "tree_covd")) {
    sim_set_sens$tau <- qweibull(0.9999, 2, 10) + 1
    sim_set_sens$model <- "tree"
  } else {
    sim_set_sens$tau <- qweibull(0.9999, 2, exp(-(-log(10) + 0 + 1/6))) + 1
    sim_set_sens$model <- "nonlinear"
  }
  for(mu in c(0.1, 0.2, 0.5)) {
    sim_set_sens$rho <- sim_set_sens$tau * mu
    simulate_LBRC_tree_methods("sensitivity analysis", sim_set_sens, current_dir)
  }
}
# Output: Supplementary Figure S13
results_dir <- file.path(current_dir, "results_intermediate")
generate_figures_tables("figure_S12_sensitivity_analysis_prediction", results_dir)


cat("\n====================================================================\n")
cat("All specified simulations and figure generations are complete.\n")
cat("Please check the 'results' directory for the output.\n")
cat("====================================================================\n")
