rm(list=ls())

source("./simulation/simulations.R")

args <- commandArgs(trailingOnly = TRUE)

M        <- as.integer(args[1])
mode_num <- args[2]

mode <- switch(as.character(mode_num),
               "1" = "test unbiasedness",
               "2" = "ANOVA study",
               "3" = "OOB tuning validation",
               "4" = "model prediction",
               "5" = "sensitivity analysis")

if(mode == "test unbiasedness"){
  Dist <- args[3]

  sim_set_list <- list(
    ksi = 500,
    M = M,
    Dist = Dist,
    ns = 200L,
    cens_rates = 20
  )
}else if(mode == "ANOVA study") {
  model <- args[3]
  Dist  <- args[4]
  n     <- as.integer(args[5])
  cens  <- as.integer(args[6])

  sim_set_list <- list(
    cov_set_num = 10,
    ksi = 500,
    M = M,
    model = model,
    Dist = Dist,
    ns = n,
    cens_rates = cens
  )
}else if(mode == "OOB tuning validation"){
  model <- args[3]
  Dist  <- args[4]
  n     <- as.integer(args[5])
  cens  <- 20
  tune.metric <- args[6]

  sim_set_list <- list(
    cov_set_num = 10,
    ksi = 500,
    M = M,
    model = model,
    Dist = Dist,
    ns = n,
    cens_rates = cens,
    tune.metric = tune.metric
  )
}else if(mode == "model prediction") {
  model <- args[3]
  Dist  <- args[4]
  n     <- as.integer(args[5])
  cens  <- as.integer(args[6])

  sim_set_list <- list(
    cov_set_num = 10,
    ksi = 500,
    M = M,
    model = model,
    Dist = Dist,
    ns = n,
    cens_rates = cens
  )
}else if(mode == "sensitivity analysis") {
  scenario <- args[3]
  mu       <- as.numeric(args[4])

  sim_set_list <- list(
    ns = 200L,
    cens_rates = 20,
    Dist = "WI",
    scenario = scenario,
    M = M
  )

  if (scenario %in% c("unbias_texpt", "unbias_covd")) {
    sim_set_list$tau <- qweibull(0.9999, 2, 3) + 1
  } else if (scenario %in% c("tree_texpt", "tree_covd")) {
    sim_set_list$cov_set_num <- 10
    sim_set_list$model <- "tree"
    sim_set_list$tau <- qweibull(0.9999, 2, 10) + 1

  } else if (scenario %in% c("nlin_texpt", "nlin_covd")) {
    sim_set_list$cov_set_num <- 10
    sim_set_list$model <- "nonlinear"
    sim_set_list$tau <- qweibull(0.9999,2,exp(-(-log(10) + 0 + 1/6))) + 1
  } else {
    stop("Unknown sensitivity scenario: ", scenario)
  }
  sim_set_list$rho <- sim_set_list$tau * mu
} else {
  stop("Unknown mode: ", mode)
}

simulate_LBRC_tree_methods(
  mode,
  sim_set_list,
  getwd()
)
