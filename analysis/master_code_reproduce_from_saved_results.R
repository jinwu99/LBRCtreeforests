rm(list=ls())

############################################################
## Master script: Reproduce all figures from saved results
##
## This script Uses generate_figs_tabs.R to regenerate ALL figures
## in the main manuscript and Supplementary Material that are
## based on simulations.
##
## IMPORTANT:
##  - This script is the recommended entry point for fast
##    reproduction of the figures without re-running any
##    simulations.
##  - For full recomputation of simulation results, use:
##      master_code_recompute_all_simulations.R
############################################################


# ==============================================================================
# Package Check & Installation
# ==============================================================================
# This script automatically checks whether all required packages
# are installed and installs any missing packages from CRAN.
install_if_missing <- function(packages) {
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
  
  if (length(to_install) > 0) {
    message("Installing missing packages: ",
            paste(to_install, collapse = ", "))
    install.packages(to_install, dependencies = TRUE)
  } else {
    message("All required packages are already installed.")
  }
}

required_packages <- c(
  "survival",
  "partykit",
  "ipred",
  "prodlim",
  "Rcpp",
  "rsample",
  "dplyr",
  "purrr",
  "tidyr",
  "ggplot2",
  "ggforce",
  "ggpubr",
  "simplecolors",
  "stringr",
  "emmeans",
  "patchwork",
  "broom",
  "scales",
  "magick",
  "pdftools"
)

install_if_missing(required_packages)


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
source("./simulation/simulations.R")
source("./simulation/generate_figs_tabs.R")
source("./simulation/real_data_application.R")
results_dir <- paste0(current_dir, "/results")


# Figure 2: recovery rate comparison across tree models
mode <- "figure_2_compare_recovery_rate"
generate_figures_tables(mode, results_dir)


# Figure 3: efficiency gain of LBRCtrees over LTRCtree
mode <- "figure_3_LBRCtrees_vs_LTRCtree_ANOVA"
generate_figures_tables(mode, results_dir)


# Figure 4: validation of OOB tuning for LBRCforests
mode <- "figure_4_WI_LBRCforests_OOB_tuning_brier"
generate_figures_tables(mode, results_dir)


# Figure 5: prediction comparison across models
mode <- "figure_5_WI_20_compare_prediction_accuracy"
generate_figures_tables(mode, results_dir)


# Figure 6~8, Table 1: real data application summary
mode <- "figure_678_table_1_real_data_application"
generate_figures_tables(mode, results_dir)


# Supplementary figure S1: test of unbiasedness of variable selection of LBRCtrees
mode <- "figure_S1_test_unbiasedness_LBRCtrees"
generate_figures_tables(mode, results_dir)


# Supplementary figure S2 ~ S6: Additional OOB tuning results
for (mode in c("figure_S2_WD_LBRCforests_OOB_tuning_brier",
               "figure_S3_LgnBat_LBRCforests_OOB_tuning_brier",
               "figure_S4_WI_LBRCforests_OOB_tuning_cindex",
               "figure_S5_WD_LBRCforests_OOB_tuning_cindex",
               "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")) {
  generate_figures_tables(mode, results_dir)
}


# Supplementary figures S7 ~ S10: Additional prediction accuracy plots
for (mode in c("figure_S7_WI_50_compare_prediction_accuracy",
               "figure_S8_WD_20_compare_prediction_accuracy",
               "figure_S9_WD_50_compare_prediction_accuracy",
               "figure_S10_LgnBat_2050_compare_prediction_accuracy")) {
  generate_figures_tables(mode, results_dir)
}


# Supplementary figure S11: sensitivity analysis on unbiasedness of variable selection
mode <- "figure_S11_sensitivity_analysis_unbiasedness"
generate_figures_tables(mode, results_dir)


# Supplementary figures S12: sensitivity analysis on tree recovery and prediction
mode <- "figure_S12_sensitivity_analysis_prediction"
generate_figures_tables(mode, results_dir)







