## Repository Structure & Folder Overview

### Root Level Scripts
Located in the same directory as this README, these two master scripts serve as the primary entry points to run the entire project or instantly reproduce the results.
* `master_code_reproduce_from_saved_results.R`: Runs the fast reproduction workflow using pre-computed outputs (detailed in the next section).
* `master_code_recompute_all_simulations.R`: Executes the full or modular recomputation workflow from scratch (detailed in the next section).

### Subdirectories

#### 1. data

Provides data-generation utilities for simulation studies, including length-biased right-censored (LBRC) sampling routines and functions for constructing various underlying survival models.

#### 2. methods

Contains the methodology implementation for LBRC conditional inference trees and forests (LBRC-CIT/CIF), including model fitting, hyperparameter tuning, prediction, performance metrics, and visualization utilities. Most components are adapted and extended from: https://github.com/weichiyao/TimeVaryingData_LTRCforests/tree/main

#### 3. metric

Provides performance evaluation functions for simulation studies, including integrated $L^2$​ prediction error and tree-structure recovery metrics.

#### 4. figures

Stores all figures included in the main manuscript and Supplementary Material, including real-data analysis outputs.

#### 5. results
Stores all pre-computed simulation results and compiled data summaries used exclusively to instantly regenerate manuscript figures without running new simulations (this fast reproducibility workflow is detailed in the next section). 

The folder structure is organized as follows:
* `properties_LBRC-CITs/`: Contains pre-computed outputs corresponding to Section 3.1 of the main manuscript and Section B of the Supplementary Material.
* `properties_LBRC-CIFs/`: Contains pre-computed outputs corresponding to Section 3.2 of the main manuscript and Sections C and D of the Supplementary Material.
* `sensitivity_analysis/`: Contains pre-computed outputs corresponding to Section E of the Supplementary Material.

Each subfolder contains:
* Pre-computed outputs organized by the underlying regression structure (`tree/`, `linear/`, `nonlinear/`, `interaction/`), with further subdirectories specified by failure time distributions and sample sizes (e.g., `tree/WI/N200/`).
* Individual `.RData` files following the standard naming convention: `LBRC_DIST_<DIST>_MODEL_<MODEL>_P<P>_N<N>_C<C>`, which explicitly indicates the failure distribution, model structure, number of covariates, sample size, and censoring rate.
* `test_unbiasedness/`: Contains pre-computed simulation results dedicated to validating the unbiased variable selection properties of LBRC-CIT estimators.

#### 6. simulation
Contains the main functions for running the simulation studies and regenerating all manuscript figures/tables. These scripts execute LBRC-CIT/CIF model training across different data-generating scenarios and construct the processed result summaries used for final figure generation.
* *Note on Intermediate Results:* When executing the full recomputation scripts (detailed below), a `results_intermediate/` directory will be automatically generated to store all newly computed `.RData` files.



## Reproducibility of Simulation Studies

We provide two distinct workflows to reproduce the simulation results, depending on your available time constraints and computational resources.

| Option | Approach                                   | Target Script                                |
| :----- | :----------------------------------------- | :------------------------------------------- |
| **A**  | **Fast Reproduction** (From saved results) | `master_code_reproduce_from_saved_results.R` |
| **B**  | **Modular Recomputation** (From scratch)   | `master_code_recompute_all_simulations.R`    |

### Option A. Fast Reproduction

Regenerates all figures and tables in the main manuscript and Supplementary Material using pre-computed intermediate results, without running any new simulations.

* **Instruction**: Open and run `master_code_reproduce_from_saved_results.R` in your R IDE. Figures will be automatically saved under the `results/` directory.

### Option B. Modular Recomputation

To accommodate reviewers' time constraints and specific areas of interest, the full recomputation script (`master_code_recompute_all_simulations.R`) is modularized. This script allows you to either execute the entire workflow sequentially or selectively run specific sections to generate target figures.

The code blocks outlined in the sections below are provided to give a transparent, step-by-step overview of the underlying simulation logic. For actual computation, we recommend opening the master script and executing the corresponding self-contained sections directly.

---

#### 1. Setup & User Settings (⚠️ MUST RUN FIRST)

Before running any specific section, you must run this setup block to initialize the environment, load required functions, and set the simulation size. 

> **Tip for Fast Reproducibility:** By default, `M_pred` and `M_test` are set to large numbers (500 and 10000), which may take several days. The code below temporarily reduces these parameters for a quick verification.

```R
# 0. Setup Directory and Functions
get_current_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) return(dirname(rstudioapi::getActiveDocumentContext()$path))
  if (!is.null(sys.frame(1)$ofile)) return(dirname(normalizePath(sys.frame(1)$ofile)))
  return(getwd())
}
current_dir <- get_current_dir()
results_dir <- file.path(current_dir, "results")
cat("Script directory:", current_dir, "\n")
setwd(current_dir)

# Load simulation and plotting functions
source("./simulation/simulations.R")
source("./simulation/generate_figs_tabs.R")

# 1. USER SETTINGS (Adjustable for fast reproduction)
M_pred <- 30    # Reduced from 500
M_test <- 1000   # Reduced from 10000
```

#### 2. Execute by Section: Main Manuscript

After running the setup block above, you can run any of the self-contained blocks below. The target figure will be automatically saved to the `results/` directory.

**Section 3.1.2: Recovering the correct tree structure (Manuscript Figure 2)**

```R
cat("\n--- Running Section 3.1.2 ---\n")
sim_set_recovery <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                         cens_rates = c(20, 50), ksi = 500, model = "tree")
for(Dist in c("WI", "WD", "Lgn", "Bat")) {
  sim_set_recovery$Dist <- Dist
  simulate_LBRC_tree_methods("ANOVA study", sim_set_recovery, current_dir)
}
generate_figures_tables("figure_2_compare_recovery_rate", results_dir)
```

**Section 3.1.3: Prediction accuracy against LTRC-CIT (Manuscript Figure 3)**

```R
cat("\n--- Running Section 3.1.3 ---\n")
sim_set_anova <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                      cens_rates = c(20, 50), ksi = 500)
for(Dist in c("WI", "WD", "Lgn", "Bat")) {
  sim_set_anova$Dist <- Dist
  if(Dist %in% c("WI", "WD")){
    for(model in c("tree", "linear", "nonlinear", "interaction")) {
      sim_set_anova$model <- model
      simulate_LBRC_tree_methods("ANOVA study", sim_set_anova, current_dir)
    }
  }else{
    sim_set_anova$model <- "tree"
    simulate_LBRC_tree_methods("ANOVA study", sim_set_anova, current_dir)
  }
}
generate_figures_tables("figure_3_LBRCtrees_vs_LTRCtree_ANOVA", results_dir)
```

**Section 3.2.1: Regulating the construction of trees in forests - WI Setting (Manuscript Figure 4)**

```R
cat("\n--- Running Section 3.2.1 ---\n")
sim_set_tune_main <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                          cens_rates = 20, ksi = 500, tune.metric = "brier")
for(model in c("tree", "linear", "nonlinear", "interaction")) {
  sim_set_tune_main$model <- model
  sim_set_tune_main$Dist <- "WI"
  simulate_LBRC_tree_methods("OOB tuning validation", sim_set_tune_main, current_dir)
}
generate_figures_tables("figure_4_WI_LBRCforests_OOB_tuning_brier", results_dir)
```

**Section 3.2.2: Prediction accuracy across methods - WI Setting, 20% Censoring (Manuscript Figure 5)**

```R
cat("\n--- Running Section 3.2.2 ---\n")
sim_set_pred_main <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400),
                          cens_rates = 20, ksi = 500)
for(model in c("tree", "linear", "nonlinear", "interaction")) {
  sim_set_pred_main$model <- model
  sim_set_pred_main$Dist <- "WI"
  simulate_LBRC_tree_methods("model prediction", sim_set_pred_main, current_dir)
}
generate_figures_tables("figure_5_WI_20_compare_prediction_accuracy", results_dir)
```

#### 3. Execute by Section: Supplementary Material

**Section B: Test of unbiasedness of variable selection (Supp Figure S1)**

```R
cat("\n--- Running Section B ---\n")
sim_set_unbias <- list(M = M_test, ns = 200L, cens_rates = c(20, 50), ksi = 500)
for(Dist in c("WI", "WD", "Lgn")) {
  sim_set_unbias$Dist <- Dist
  simulate_LBRC_tree_methods("test unbiasedness", sim_set_unbias, current_dir)
}
generate_figures_tables("figure_S1_test_unbiasedness_LBRCtrees", results_dir)
```

**Section C.2.1: Results based on IBS (Supp Figures S2-S3)**

```R
cat("\n--- Running Section C.2.1 ---\n")
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
for (fig in c("figure_S2_WD_LBRCforests_OOB_tuning_brier", "figure_S3_LgnBat_LBRCforests_OOB_tuning_brier")) {
  generate_figures_tables(fig, results_dir)
}
```

**Section C.2.2: Results based on C-index (Supp Figures S4-S6)**

```R
cat("\n--- Running Section C.2.2 ---\n")
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
for (fig in c("figure_S4_WI_LBRCforests_OOB_tuning_cindex", "figure_S5_WD_LBRCforests_OOB_tuning_cindex", "figure_S6_LgnBat_LBRCforests_OOB_tuning_cindex")) {
  generate_figures_tables(fig, results_dir)
}
```

**Section D: Additional prediction accuracy across methods (Supp Figures S7-S10)**

```R
cat("\n--- Running Section D ---\n")
sim_set_pred_supp <- list(M = M_pred, cov_set_num = 10, ns = c(100, 200, 400), ksi = 500)
for(dist in c("WI", "WD", "Lgn", "Bat")){
  if(dist == "WI"){
    for(model in c("tree", "linear", "nonlinear", "interaction")) {
      sim_set_pred_supp$model <- model; sim_set_pred_supp$Dist <- dist; sim_set_pred_supp$cens_rates <- 50
      simulate_LBRC_tree_methods("model prediction", sim_set_pred_supp, current_dir)
    }
  }else if(dist == "WD"){
    for(model in c("tree", "linear", "nonlinear", "interaction")) {
      sim_set_pred_supp$model <- model; sim_set_pred_supp$Dist <- dist; sim_set_pred_supp$cens_rates <- c(20, 50)
      simulate_LBRC_tree_methods("model prediction", sim_set_pred_supp, current_dir)
    }
  }else{
    sim_set_pred_supp$model <- "tree"; sim_set_pred_supp$Dist <- dist; sim_set_pred_supp$cens_rates <- c(20, 50)
    simulate_LBRC_tree_methods("model prediction", sim_set_pred_supp, current_dir)
  }
}
for (fig in c("figure_S7_WI_50_compare_prediction_accuracy", "figure_S8_WD_20_compare_prediction_accuracy", "figure_S9_WD_50_compare_prediction_accuracy", "figure_S10_LgnBat_2050_compare_prediction_accuracy")) {
  generate_figures_tables(fig, results_dir)
}
```

**Section E.1: Unbiasedness of variable selection (Sensitivity) (Supp Figure S12)**

```R
cat("\n--- Running Section E.1 ---\n")
sim_set_sens <- list(ns = 200L, cens_rates = 20, Dist = "WI")
for(scenario in c("unbias_texpt", "unbias_covd")) {
  sim_set_sens$scenario <- scenario; sim_set_sens$M <- M_test; sim_set_sens$tau <- qweibull(0.9999, 2, 3) + 1
  for(mu in c(0.1, 0.2, 0.5)) {
    sim_set_sens$rho <- sim_set_sens$tau * mu
    simulate_LBRC_tree_methods("sensitivity analysis", sim_set_sens, current_dir)
  }
}
generate_figures_tables("figure_S12_sensitivity_analysis_unbiasedness", results_dir)
```

**Section E.2: Tree recovery and prediction accuracy (Sensitivity) (Supp Figure S13)**

```R
cat("\n--- Running Section E.2 ---\n")
sim_set_sens <- list(ns = 200L, cens_rates = 20, Dist = "WI")
for(scenario in c("tree_texpt", "tree_covd", "nlin_texpt", "nlin_covd")) {
  sim_set_sens$scenario <- scenario; sim_set_sens$M <- M_pred; sim_set_sens$cov_set_num <- 10
  if(scenario %in% c("tree_texpt", "tree_covd")) {
    sim_set_sens$tau <- qweibull(0.9999, 2, 10) + 1; sim_set_sens$model <- "tree"
  } else {
    sim_set_sens$tau <- qweibull(0.9999, 2, exp(-(-log(10) + 0 + 1/6))) + 1; sim_set_sens$model <- "nonlinear"
  }
  for(mu in c(0.1, 0.2, 0.5)) {
    sim_set_sens$rho <- sim_set_sens$tau * mu
    simulate_LBRC_tree_methods("sensitivity analysis", sim_set_sens, current_dir)
  }
}
generate_figures_tables("figure_S13_sensitivity_analysis_prediction", results_dir)
```

