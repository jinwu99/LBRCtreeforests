## Reproducibility Instructions
Two master scripts are provided for reproducing simulation-based results:
### 1. Fast reproduction of all manuscript figures
* `master_code_reproduce_from_saved_results.R` 
  Regenerates all figures in the main manuscript and Supplementary Material using pre-computed results, without running any simulations. Figures are automatically saved under the `results/` directory.
### 2. Full recomputation of simulations
* `master_code_reproduce_from_saved_results.R`
  Re-runs **all** simulations used in the paper and regenerates all intermediate `.RData` files. A `results_intermediate/` directory is automatically created to store newly generated simulation outputs.
  * Warning: Full recomputation requires **multiple days** of runtime on a standard computer. For quick reproducibility checks, please use the fast script above.

## Folder Overview
### 1. data
Provides data-generation utilities for simulation studies, including length-biased right-censored (LBRC) sampling routines and functions for constructing various underlying survival models.
### 2. methods
Contains the methodology implementation for LBRC conditional inference trees and forests (LBRC-CIT/CIF), including model fitting, hyperparameter tuning, prediction, performance metrics, and visualization utilities. Most components are adapted and extended from: https://github.com/weichiyao/TimeVaryingData_LTRCforests/tree/main
### 3. metric
Provides performance evaluation functions for simulation studies, including integrated $L^2$ prediction error and tree-structure recovery metrics.
### 4. results
Stores all pre-computed simulation results and manuscript figures used for fast reproducibility. This includes:
- Simulation result files organized by data-generating setting (`tree/`, `linear/`, `nonlinear/`), with subfolders by distribution and sample size (e.g., `tree/WI/N200/`). Each `.RData` file follows the naming convention:  `LBRC_DIST_<DIST>_MODEL_<MODEL>_P<P>_N<N>_C<C>` indicating failure distribution, model structure, number of covariates, sample size, and censoring rate.
- `test_unbiasedness/` containing simulation results for validating unbiasedness of LBRC-CIT estimators.
- `figures/` containing all figures and tables included in the main manuscript and Supplementary Material, including real-data analysis outputs.
### 5. simulation
Contains main functions for running the simulation studies and regenerating all manuscript figures/tables. These scripts execute LBRC-CIT/CIF model training across different data-generating scenarios, save results to the `results_intermediate/` directory, and construct the processed result summaries used for figure and table generation.
