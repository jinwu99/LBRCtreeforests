rm(list=ls())

############################################################
## Master script: Reproduce all figures from saved results
##
## This script Uses plot_results.R to regenerate ALL figures
## in the main manuscript and Supplementary Material that are
## based on simulations.
##
## IMPORTANT:
##  - This script is the recommended entry point for fast
##    reproduction of the figures without re-running any
##    simulations.
##  - For full recomputation of simulation results, use:
##      master_recompute_all_simulations.R
############################################################

## Utility: set working directory to the location of this script ----
get_current_dir <- function() {
  # Case 1: running in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }

  # Case 2: running via source() or in R console
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }

  # Case 3: running as an Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }

  # Fallback: working directory
  return(getwd())
}


current_dir <- get_current_dir()
cat("Script directory:", current_dir, "\n")
setwd(current_dir)
source("./simulation/generate_figs_tabs.R")
results_dir <- paste0(current_dir, "/results")


# Figure 1: recovery rate comparison across tree models
mode <- "figure_1_compare_recovery_rate"
generate_figures_tables(mode, results_dir)


# Figure 2: efficiency gain of LBRCtrees over LTRCtree
mode <- "figure_2_LBRCtrees_vs_LTRCtree_ANOVA"
generate_figures_tables(mode, results_dir)


# Figure 3: validation of OOB tuning for LBRCforests
mode <- "figure_3_LBRCforests_OOB_tuning"
generate_figures_tables(mode, results_dir)


# Figure 4: prediction comparison across models
mode <- "figure_4_compare_prediction_accuracy"
generate_figures_tables(mode, results_dir)


# Supplementary table S1: test of unbiasedness of variable selection of LBRCtrees
mode <- "table_S1_test_unbiasedness_LBRCtrees"
generate_figures_tables(mode, results_dir)


# Supplementary figure S1: Additional OOB tuning results
mode <- "figure_S1_LBRCforests_OOB_tuning"
generate_figures_tables(mode, results_dir)


# Supplementary figures S2 and S3: Additional prediction accuracy plots
for (mode in c("figure_S2_compare_prediction_accuracy",
               "figure_S3_compare_prediction_accuracy")) {
  generate_figures_tables(mode, results_dir)
}







