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



# ------------ Main simulation: CIT / CIF performance studies ------------
# These settings correspond to Section 3 and Supplements:
# - Tree structure recovery
# - Prediction accuracy of LBRC-CIT vs LTRC-CIT (ANOVA-type targeting unbiased survival function)
# - CIF tuning and prediction comparisons (targeting observable conditional survival function)
# - Sensitivity analysis to the violation of stationary assumption on truncation process


# Simulation modes:
# - "ANOVA stuudy":           tree
# - "model prediction":       proposed CIT/CIF prediction accuracy
# - "OOB tuning validation":  mtry tuning via OOB IBS
# - "test unbiasedness":      check unbiased variable selection in CIT
# - "sensitivity analysis":   sensitivity analysis
simulation_mode <- c("ANOVA study",
                     "model prediction",
                     "OOB tuning validation",
                     "test unbiasedness",
                     "sensitivity analysis")

# Total number of simulations for "model prediction" and "OOB tuning validation"
M_pred <- 500

# Total number of simulations for "test unbiasedness"
M_test <- 10000


# Loop over simulation mode, underlying regression structure, and
# failure time distribution (see Section 3 and figure captions).
# - WI:  Weibull with increasing hazard
# - WD:  Weibull with decreasing hazard
# - Lgn: Lognormal
# - Bat: Bathtub-shaped hazard
for(mode in simulation_mode) {
  if(mode %in% c("ANOVA study", "model prediction", "OOB tuning validation")){
    # Simulation design list
    sim_set_list <- list()
    # Total number of simulations
    sim_set_list$M <- M_pred
    # Number of covariate sets; each set has three types (continuous, ordinal, binary)
    sim_set_list$cov_set_num <- 10
    # Sample sizes considered in the paper
    sim_set_list$ns <- c(100, 200, 400)

    if(mode == "OOB tuning validation"){
      # Censoring rates (in percentages)
      sim_set_list$cens_rates <- 20
      # tuning metric
      sim_set_list$tune.metric <- "cindex"
    }else{
      sim_set_list$cens_rates <- c(20, 50)
    }

    # Upper bound for truncation time; set large for length-biased sampling assumption
    sim_set_list$ksi <- 500

    for(model in c("tree", "linear", "nonlinear", "interaction")){
      sim_set_list$model <- model
      if(model == "tree"){
        for(Dist in c("WI", "WD", "Lgn", "Bat")){
          sim_set_list$Dist <- Dist
          simulate_LBRC_tree_methods(mode, sim_set_list, current_dir)
        }
      }else{
        for(Dist in c("WI", "WD")){
          sim_set_list$Dist <- Dist
          simulate_LBRC_tree_methods(mode, sim_set_list, current_dir)
        }
      }
    }
  }else if(mode == "test unbiasedness"){
    sim_set_list <- list()
    sim_set_list$M <- M_test
    sim_set_list$ns <- 200L
    sim_set_list$cens_rates <- c(20, 50)
    sim_set_list$ksi <- 500

    for(Dist in c("WI", "WD", "Lgn")){
      sim_set_list$Dist <- Dist
      simulate_LBRC_tree_methods(mode, sim_set_list, current_dir)
    }
  }else if(mode == "sensitivity analysis"){
    sim_set_list <- list()
    sim_set_list$ns <- 200L
    sim_set_list$cens_rates <- 20
    sim_set_list$Dist <- "WI"

    for(scenario in c("unbias_texpt", "unbias_covd", "tree_texpt", "tree_covd", "nlin_texpt", "nlin_covd")){
      sim_set_list$scenario <- scenario
      if(scenario %in% c("unbias_texpt", "unbias_covd")){
        sim_set_list$M <- M_test
        sim_set_list$tau <- qweibull(0.9999,2,3) + 1
      }else{ # scenario %in% c("tree_texpt", "tree_covd", "nlin_texpt", "nlin_covd")
        sim_set_list$M <- M_pred
        sim_set_list$cov_set_num <- 10
        if(scenario %in% c("tree_texpt", "tree_covd")){
          sim_set_list$tau <- qweibull(0.9999,2,10) + 1
          sim_set_list$model <- "tree"
        }else{ # scenario %in% c("nlin_texpt", "nlin_covd")
          sim_set_list$tau <- qweibull(0.9999,2,exp(-(-log(10) + 0 + 1/6))) + 1
          sim_set_list$model <- "nonlinear"
        }
      }
      for(mu in c(0.1, 0.2 ,0.5)){ # mu = rho / tau
        sim_set_list$rho <- sim_set_list$tau * mu
        simulate_LBRC_tree_methods(mode, sim_set_list, current_dir)
      }
    }
  }
}



##########################################################################
## Plot all the results ----------

current_dir <- get_current_dir()
cat("Script directory:", current_dir, "\n")
setwd(current_dir)
source("./simulation/generate_figs_tabs.R")
results_dir <- paste0(current_dir, "/results")


## Figure 1: recovery rate comparison across tree models
mode <- "figure_1_compare_recovery_rate"
generate_figures_tables(mode, results_dir)


# Figure 2: efficiency gain of LBRCtrees over LTRCtree
mode <- "figure_2_LBRCtrees_vs_LTRCtree_ANOVA"
generate_figures_tables(mode, results_dir)


# Figure 3: validation of OOB tuning for LBRCforests
mode <- "figure_3_WI_LBRCforests_OOB_tuning_brier"
generate_figures_tables(mode, results_dir)


# Figure 4: prediction comparison across models
mode <- "figure_4_WI_20_compare_prediction_accuracy"
generate_figures_tables(mode, results_dir)



# Supplementary figure S1: test of unbiasedness of variable selection of LBRCtrees
mode <- "figure_S1_test_unbiasedness_LBRCtrees"
generate_figures_tables(mode, results_dir)


# Supplementary figure S2 ~ S6: Additional OOB tuning results
for (mode in c("figure_S2_WI_LBRCforests_OOB_tuning_cindex",
               "figure_S3_WD_LBRCforests_OOB_tuning_brier",
               "figure_S4_WD_LBRCforests_OOB_tuning_cindex",
               "figure_S5_LgnBat_LBRCforests_OOB_tuning_brier",
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





