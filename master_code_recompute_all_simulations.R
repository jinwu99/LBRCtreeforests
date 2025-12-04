rm(list=ls())

############################################################
## Master script: Recompute all simulation results
##
## This script reruns ALL simulation studies from scratch
## and regenerates the intermediate results.
##
## IMPORTANT:
##  - This script is intended for full recomputation and may
##    take several days to finish.
##  - For fast reproduction of figures/tables from saved results,
##    use: "master_code_reproduce_from_saved_results.R".
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
source("./simulation/simulations.R")


## Main simulation: CIT / CIF performance studies ------------
# These settings correspond to Section 3.1 and 3.2 of the paper:
# - Tree structure recovery
# - Prediction accuracy of LBRC-CIT vs LTRC-CIT (ANOVA-type)
# - CIF tuning and prediction comparisons


# Simulation modes:
# - "model prediction":       tree recovery + prediction accuracy
# - "OOB tuning validation":  mtry tuning via OOB IBS
# - "test unbiasedness":      check unbiased variable selection in CIT
simulation_mode <- c("model prediction", "OOB tuning validation", "test unbiasedness")

# Total number of simulations for "model prediction" and "OOB tuning validation"
M_pred <- 500L

# Total number of simulations for "test unbiasedness"
M_test <- 10000L


# Loop over simulation mode, underlying regression structure, and
# failure time distribution (see Section 3 and figure captions).
# - WI:  Weibull with increasing hazard
# - WD:  Weibull with decreasing hazard
# - Lgn: Lognormal
# - Bat: Bathtub-shaped hazard
for(mode in simulation_mode) {
  if(mode %in% c("model prediction", "OOB tuning validation")){
    # Total number of simulations
    M <- M_pred
    # Simulation design list
    sim_set_list <- list()
    # Number of covariate sets; each set has three types (continuous, ordinal, binary)
    sim_set_list$cov_set_num <- 10L
    # Sample sizes considered in the paper
    sim_set_list$ns <- c(100L, 200L, 400L)
    # Censoring rates (in percentages)
    sim_set_list$cens_rates <- c(20, 50)
    # Upper bound for truncation time; set large for length-biased sampling assumption
    sim_set_list$ksi <- 500

    for(model in c("tree", "linear", "nonlinear")){
      sim_set_list$model <- model
      if(model == "tree"){
        for(Dist in c("WI", "WD", "Lgn", "Bat")){
          sim_set_list$Dist <- Dist
          simulate_LBRC_tree_methods(mode, sim_set_list, M, current_dir)
        }
      }else{
        for(Dist in c("WI", "WD")){
          sim_set_list$Dist <- Dist
          simulate_LBRC_tree_methods(mode, sim_set_list, M, current_dir)
        }
      }
    }
  }else{ # test unbiasedness
    M <- M_test
    sim_set_list <- list()
    sim_set_list$ns         <- 300L
    sim_set_list$cens_rates <- c(20, 50)
    sim_set_list$ksi        <- 100

    for(Dist in c("WI", "WD", "Lgn")){
      sim_set_list$Dist <- Dist
      simulate_LBRC_tree_methods(mode, sim_set_list, M, current_dir)
    }
  }
}


##########################################################################
## Plot all the results ----------

current_dir <- get_current_dir()
cat("Script directory:", current_dir, "\n")
setwd(current_dir)
source("./simulation/generate_figs_tabs.R")
results_dir <- paste0(current_dir, "/results_intermediate")


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


