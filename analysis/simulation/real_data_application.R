###############################################################
## Real-data application
##
## Code for the real-data analyses
##   - stationarity assessment,
##   - fitted survival trees,
##   - repeated cross-validation,
##   - Brier score trajectories,
##   - integrated Brier score summaries.
##
###############################################################

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(magick)
library(pdftools)

# Formal test of the stationarity assumption in prevalent cohort studies.
# Adapted from the archived CoxPhLb package:
# https://cran.r-project.org/src/contrib/Archive/CoxPhLb/
# Reference: Addona and Wolfson (2006).
station.test <- function(a, v, delta, digits = 3L)
{
  len<-length(delta)

  sum<-0

  for(i in 1:len)
    for(j in 1:len)
    {
      sum<- sum + ifelse(a[i]>v[j] & delta[j]==1, 1, 0) - ifelse(a[i]<v[j], 1, 0)
    }
  # Calculate test statistic

  wtest <- sum/len^2
  # Calculate the variance of the test statistic
  var <- 0
  for( j in 1:len)
  {
    var <- var + sum(ifelse(v>a[j], 1, 0))^2+
      sum(ifelse(a>v[j]& delta[j]==1, 1, 0))^2+
      2*sum(ifelse(v<=a[j] & delta==1, 1, 0))*sum(ifelse(a>v[j], 1, 0))*delta[j]+
      2*sum(ifelse(v>a[j], 1, 0))*sum(ifelse(a<=v[j], 1, 0))+
      2*sum(ifelse(a<=v[j], 1, 0))*sum(ifelse(v<=a[j] & delta==1, 1, 0))-
      2*sum(ifelse(a>v[j], 1, 0))*sum(ifelse(v>a[j], 1, 0))*delta[j]

  }
  # standardized test statistics
  test.statistic <- sqrt(len)*wtest/sqrt(var/len^3)
  p.value <- 1 + pnorm(-abs(test.statistic)) - pnorm(abs(test.statistic))

  tab = cbind(test.statistic = format(round(test.statistic, digits)), p.value = format.pval(p.value, digits, eps = 1e-03))
  rownames(tab) = " "
  tab = as.table(tab)
  print(tab)
  out=list(test.statistic = test.statistic, p.value = p.value, result = tab)
  class(out) <- "station.test"
  invisible(out)
}
station.test.plot <- function(a, v, delta)
{
  len<-length(delta)

  plot(survfit(Surv(v,delta)~1), conf.int = FALSE, lty=1, ylab="Estimated survival function")
  lines(survfit(Surv(a,rep(1,len))~1), conf.int = FALSE, lty=2)
  legend("topright", c("forward recurrence time","backward recurrence time"), lty=c(1,2))
}

# Compute risk-set sizes over follow-up time and identify
# the last time point at which a specified percentage of
# subjects remain at risk.
get_risk_set_info <- function(DATA, target_pct) {
  N_total <- nrow(DATA)
  unique_event_times <- sort(unique(DATA$Y[DATA$event == 1]))
  risk_set_df <- tibble(Time = unique_event_times) %>%
    rowwise() %>%
    mutate(N_Risk = sum(DATA$Y >= Time), Pct_Risk = (N_Risk / N_total) * 100) %>%
    ungroup()
  
  valid_times <- risk_set_df %>% filter(Pct_Risk >= target_pct)
  target_time_point <- ifelse(nrow(valid_times) > 0, max(valid_times$Time), NA)
  
  return(list(risk_set_df = risk_set_df, target_time = target_time_point))
}

# Real-data prediction analysis using repeated cross-validation.
# Computes Brier score curves and integrated Brier scores (IBS)
# for all competing LTRC and LBRC methods and reproduces
# the prediction results reported in the manuscript.
cross_validation_prediction <- function(DATA,
                                        covariates,
                                        v = 10,
                                        repeats = 10,
                                        ntree = 100,
                                        batch_size = 50,
                                        time.tau = NULL,
                                        out_dir = ".") {
  synth_data_location <- paste0(out_dir, "/synth_data_application")
  if (!dir.exists(synth_data_location)) {
    dir.create(synth_data_location, recursive = TRUE)
  }

  DATA$id <- seq_len(nrow(DATA))

  formula <- as.formula(
    paste("Surv(A, Y, event) ~", paste(covariates, collapse = " + "))
  )

  n <- nrow(DATA)

  vardi_tune <- list(eps = 1e-7, max_iter = 20)

  folds <- rsample::vfold_cv(DATA, v = v, repeats = repeats)

  methods <- c(
    "LTRC-CIT",
    "LBRC-CIT-C",
    "LBRC-CIT-F",
    "LTRC-CIF",
    "LBRC-CIF-C",
    "LBRC-CIF-F",
    "LTRC-COX",
    "LBRC-COX"
  )

  plot_grid <- sort(unique(DATA$Y))
  plot_grid <- plot_grid[plot_grid <= time.tau]

  ibs_path <- file.path(synth_data_location, "synth_data_ibs_results.rds")
  bs_path  <- file.path(synth_data_location, "synth_data_bs_results.rds")

  if (file.exists(ibs_path)) {
    ibs_df_saved <- readRDS(ibs_path)
    all_ibs <- split(ibs_df_saved, seq_len(nrow(ibs_df_saved)))
    cat("Loaded existing IBS results:", nrow(ibs_df_saved), "rows\n")
  } else {
    all_ibs <- list()
  }

  if (file.exists(bs_path)) {
    bs_df_saved <- readRDS(bs_path)
    all_bs <- split(bs_df_saved, seq_len(nrow(bs_df_saved)))
    cat("Loaded existing BS results:", nrow(bs_df_saved), "rows\n")
  } else {
    all_bs <- list()
  }

  if (length(all_ibs) > 0) {
    completed <- bind_rows(all_ibs) %>%
      distinct(fold_id, Method)
  } else {
    completed <- data.frame(
      fold_id = integer(),
      Method = character(),
      stringsAsFactors = FALSE
    )
  }

  for (fold_id in seq_len(nrow(folds))) {

    cat("\n==============================\n")
    cat("Fold", fold_id, "of", nrow(folds), "\n")
    cat("==============================\n")

    split <- folds$splits[[fold_id]]
    TRAIN <- rsample::analysis(split)
    TEST  <- rsample::assessment(split)

    TEST$id <- seq_len(nrow(TEST))
    # time.uniq <- sort(unique(c(0, plot_grid)))
    time.uniq <- sort(unique(c(0, TEST$Y[TEST$Y <= time.tau])))

    control <- partykit::ctree_control(
      teststat  = "quad",
      testtype  = "Univ",
      minsplit  = max(ceiling(sqrt(nrow(TRAIN))), 20),
      minbucket = max(ceiling(sqrt(nrow(TRAIN))), 7),
      saveinfo  = FALSE
    )

    for (method in methods) {

      if (nrow(completed) > 0 &&
          any(completed$fold_id == fold_id & completed$Method == method)) {
        cat("Skipping completed:", method, "fold", fold_id, "\n")
        next
      }

      cat("Running:", method, "\n")

      res <- tryCatch({

        pred <- fit_predict_one_method(
          method = method,
          formula = formula,
          covariates = covariates,
          TRAIN = TRAIN,
          TEST = TEST,
          time.uniq = time.uniq,
          time.tau = time.tau,
          ntree = ntree,
          control = control,
          vardi_tune = vardi_tune,
          batch_size = batch_size
        )

        brier <- sbrier_lbrc(
          obj  = pred$survival.obj,
          id   = pred$survival.id,
          pred = pred,
          type = "BS"
        )

        ibs <- sbrier_lbrc(
          obj  = pred$survival.obj,
          id   = pred$survival.id,
          pred = pred,
          type = "IBS"
        )

        list(brier = brier, ibs = as.numeric(ibs))

      }, error = function(e) {

        msg <- paste(
          Sys.time(),
          "fold:", fold_id,
          "method:", method,
          "error:", conditionMessage(e)
        )

        cat("ERROR:", msg, "\n")

        write(
          msg,
          file = file.path(synth_data_location, "error_log.txt"),
          append = TRUE
        )

        NULL
      })

      if (!is.null(res)) {

        ibs_new <- data.frame(
          fold_id = fold_id,
          Method = method,
          IBS = res$ibs
        )

        bs_interp <- interpolate_brier(
          brier_df = res$brier,
          plot_grid = plot_grid
        )

        bs_new <- data.frame(
          fold_id = fold_id,
          Method = method,
          Time = plot_grid,
          BScore = bs_interp
        )

        all_ibs[[length(all_ibs) + 1]] <- ibs_new
        all_bs[[length(all_bs) + 1]] <- bs_new

        completed <- bind_rows(completed, ibs_new[, c("fold_id", "Method")]) %>%
          distinct(fold_id, Method)

        ibs_tmp <- bind_rows(all_ibs)
        bs_tmp  <- bind_rows(all_bs)

        saveRDS(ibs_tmp, ibs_path)
        saveRDS(bs_tmp, bs_path)

        cat("Intermediate results saved:", method, "fold", fold_id, "\n")
      }

      gc()
    }
  }

}

# Fit a specified survival method and generate predictions
# on the test dataset. Supports LTRC/LBRC CIT, CIF, and Cox
# model variants used in the real-data application.
fit_predict_one_method <- function(method,
                                   formula,
                                   covariates,
                                   TRAIN,
                                   TEST,
                                   time.uniq,
                                   time.tau,
                                   ntree,
                                   control,
                                   vardi_tune,
                                   batch_size) {

  if (method == "LTRC-CIT") {

    obj <- lbrccit(
      formula = formula,
      data = TRAIN,
      perm_test_est = "KM"
    )

    pred <- predictProb_LBRC(
      object = obj,
      newdata = TEST,
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST)),
      pred_surv_est = "KM",
      target = "observed"
    )

  } else if (method == "LBRC-CIT-C") {

    obj <- lbrccit(
      formula = formula,
      data = TRAIN,
      perm_test_est = "MCLE"
    )

    pred <- predictProb_LBRC(
      object = obj,
      newdata = TEST,
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST)),
      pred_surv_est = "MCLE",
      target = "observed"
    )

  } else if (method == "LBRC-CIT-F") {

    obj <- lbrccit(
      formula = formula,
      data = TRAIN,
      perm_test_est = "MFLE",
      perm_test_args = vardi_tune
    )

    pred <- predictProb_LBRC(
      object = obj,
      newdata = TEST,
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST)),
      pred_surv_est = "MFLE",
      pred_surv_args = vardi_tune,
      target = "observed"
    )

  } else if (method == "LTRC-CIF") {

    obj <- lbrccif(
      formula = formula,
      data = TRAIN,
      mtry = NULL,
      perm_test_est = "KM",
      pred_surv_est = "KM",
      ntree = ntree,
      control = control,
      trace = FALSE
    )

    pred <- batched_predictProb_LBRC(
      object = obj,
      newdata = TEST,
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST)),
      pred_surv_est = "KM",
      batch_size = batch_size
    )

  } else if (method == "LBRC-CIF-C") {

    obj <- lbrccif(
      formula = formula,
      data = TRAIN,
      mtry = NULL,
      perm_test_est = "MCLE",
      pred_surv_est = "MCLE",
      ntree = ntree,
      control = control,
      trace = FALSE
    )

    pred <- batched_predictProb_LBRC(
      object = obj,
      newdata = TEST,
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST)),
      pred_surv_est = "MCLE",
      batch_size = batch_size
    )

  } else if (method == "LBRC-CIF-F") {

    obj <- lbrccif(
      formula = formula,
      data = TRAIN,
      mtry = NULL,
      perm_test_est = "MFLE",
      perm_test_args = vardi_tune,
      pred_surv_est = "MFLE",
      pred_surv_args = vardi_tune,
      ntree = ntree,
      control = control,
      trace = FALSE
    )

    pred <- batched_predictProb_LBRC(
      object = obj,
      newdata = TEST,
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST)),
      pred_surv_est = "MFLE",
      batch_size = batch_size
    )

  } else if (method == "LTRC-COX") {

    obj <- survival::coxph(formula = formula, data = TRAIN, model=T, x=T, y=T)
    sf <- survival::survfit(obj, newdata = TEST)
    sf_list <- lapply(seq_len(nrow(TEST)), function(i) sf[i])

    pred <- predict_cox(
      pred = sf_list,
      data = TEST[, c("A", "Y", "event", "id")],
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST))
    )

  } else if (method == "LBRC-COX") {

    CPDATA <- data.frame(
      A = c(TRAIN$A, (TRAIN$Y - TRAIN$A)[TRAIN$event == 1]),
      Y = c(TRAIN$Y, TRAIN$Y[TRAIN$event == 1]),
      event = c(TRAIN$event, TRAIN$event[TRAIN$event == 1])
    )
    for (cov in covariates) {
      CPDATA[[cov]] <- c(TRAIN[[cov]], TRAIN[TRAIN$event == 1, cov])
    }
    obj <- survival::coxph(formula = formula, data = CPDATA, model=T, x=T, y=T)
    sf <- survival::survfit(obj, newdata = TEST)
    sf_list <- lapply(seq_len(nrow(TEST)), function(i) sf[i])
    pred <- predict_cox(
      pred = sf_list,
      data = TEST[, c("A", "Y", "event", "id")],
      time.eval = time.uniq,
      time.tau = rep(time.tau, nrow(TEST))
    )
  } else {
    stop("Unknown method: ", method)
  }

  ## ensure required fields for sbrier_lbrc
  pred$survival.obj <- survival::Surv(TEST$A, TEST$Y, TEST$event)
  pred$survival.id <- TEST$id
  pred$survival.times <- time.uniq
  pred$survival.tau <- rep(time.tau, nrow(TEST))

  return(pred)
}

# Convert survfit() predictions from a Cox model into the
# prediction object format required for Brier score and IBS
# evaluation.
predict_cox <- function(pred, data, time.eval, time.tau = NULL) {

  n <- nrow(data)
  if (is.null(time.tau)){
    time.tau <- rep(max(time.eval), N)
  }

  pred$survival.probs <- sapply(seq_len(n), function(i) {
    .shatfunc(
      Ni = i,
      data = data,
      pred = pred,
      tpnt = time.eval,
      tau = time.tau
    )
  })

  pred$survival.times <- time.eval
  pred$survival.tau <- time.tau
  pred$survival.obj <- survival::Surv(data$A, data$Y, data$event)
  pred$survival.id <- data$id

  return(pred)
}

# Interpolate Brier score trajectories onto a common time grid
# for averaging and visualization across cross-validation folds.
interpolate_brier <- function(brier_df, plot_grid) {

  if (is.null(brier_df) || nrow(brier_df) < 2) {
    return(rep(NA_real_, length(plot_grid)))
  }

  approx(
    x = brier_df$Time,
    y = brier_df$BScore,
    xout = plot_grid,
    rule = 2
  )$y
}


###############################################################
###############################################################
# Example workflow for real-data application code
#
# This example illustrates the main steps used in the manuscript:
#   1. assess stationarity of the length-biased sampling process,
#   2. fit LTRC/LBRC conditional inference trees,
#   3. evaluate prediction performance by repeated cross-validation.
#
# The real data obtained from the Cancer Clinical Library
# Database are not publicly available due to privacy regulations.
# Therefore, a simulated LBRC dataset is used here to demonstrate
# the analysis workflow and reproduce the code structure.
###############################################################
###############################################################

real_data_application <- function(mode        = "cross_validation",
                                  working_dir = NULL,
                                  v           = 10, # v-fold parition
                                  repeats     = 1)  # number of times to repeat the v-fold partition
{ 
  rand <- 100
  set.seed(rand)
  
  # Set directory to store files
  if (mode == "cross_validation") out_dir <- paste0(working_dir, "/results_intermediate")
  else out_dir <- working_dir
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  setwd(out_dir)
  
  synth_data_location <- paste0(out_dir, "/synth_data_application")
  if (!dir.exists(synth_data_location)) {
    dir.create(synth_data_location, recursive = TRUE)
  }
  
  # Generate an example LBRC dataset for demonstration
  synth_data_path <- file.path(paste0(out_dir,"/synth_data_application"), "synthetic_data.csv")
  if (file.exists(synth_data_path)) {
    DATA <- read.csv(synth_data_path)
  } else {
    DATA <- LBRC.generate_PH_nPH(n           = 400,
                                 Dist        = "WI",
                                 cens_rate   = 50,
                                 ksi         = 500,
                                 cov_set_num = 12,
                                 model       = "nonlinear")$Data
    write.csv(DATA, synth_data_path, row.names = FALSE)
  }
  
  # Set the prediction horizon to the last time point with at least 10% at risk
  time.tau <- get_risk_set_info(DATA = DATA, target_pct = 10)$target_time
  
  # Define covariates and survival formula
  covariates <- setdiff(colnames(DATA), c("A", "Y", "event", "id"))
  formula <- as.formula(
    paste("Surv(A, Y, event) ~", paste(covariates, collapse = " + "))
  )
  
  if(mode == "generate_figs_tabs"){
    # Check the stationarity assumption using backward and forward recurrence times
    station.test(DATA$A, (DATA$Y-DATA$A), DATA$event)
    cairo_pdf("figure_6_synth_data_stationarity.pdf", width = 8, height = 6)
    station.test.plot(DATA$A, (DATA$Y-DATA$A), DATA$event)
    dev.off()
    
    # Fit and visualize LTRC-CIT and LBRC-CIT trees
    # using the KM, MCLE, and MFLE survival estimators.
    cairo_pdf("synth_data_LTRC-CIT_result.pdf", width = 8, height = 6)
    obj <- lbrccit(
      formula = formula,
      data = DATA,
      perm_test_est = "KM"
    )
    plot(obj)
    dev.off()
    
    cairo_pdf("synth_data_LBRC-CIT-C_result.pdf", width = 8, height = 6)
    obj <- lbrccit(
      formula = formula,
      data = DATA,
      perm_test_est = "MCLE"
    )
    plot(obj)
    dev.off()
    
    cairo_pdf("synth_data_LBRC-CIT-F_result.pdf", width = 8, height = 6)
    obj <- lbrccit(
      formula = formula,
      data = DATA,
      perm_test_est = "MFLE"
    )
    plot(obj)
    dev.off()
    
    tmp_files <- c("synth_data_LTRC-CIT_result.pdf", "synth_data_LBRC-CIT-C_result.pdf", "synth_data_LBRC-CIT-F_result.pdf")
    suppressMessages({
      imgs <- lapply(tmp_files, image_read_pdf, density = 600)
    })
    imgs <- do.call(c, imgs)
    grid_3x1 <- image_append(imgs, stack = TRUE)
    image_write(grid_3x1, path = "figure_7_synth_data_CIT_results.pdf", format = "pdf")
    file.remove(tmp_files)
    cat("Successfully generated and merged tree plots into 'figure_7_synth_data_CIT_results.pdf'\n")
  }

  if(mode == "cross_validation"){
    # Run repeated cross-validation for LTRC/LBRC CIT, CIF, and Cox variants
    res <- cross_validation_prediction(
      DATA = DATA,
      covariates = covariates,
      v = v,
      repeats = repeats,
      ntree = 100,
      batch_size = 50,
      time.tau = time.tau,
      out_dir = out_dir
    )
  }
}

