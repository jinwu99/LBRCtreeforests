# Implementation adapted from LTRCforests:
# https://github.com/weichiyao/TimeVaryingData_LTRCforests/blob/main/pkg/LTRCforests/R/tune.ltrccif.R}
tune.lbrccif  <- function(formula, data, id,
                          perm_test_est = "MCLE", perm_test_args = list(),
                          pred_surv_est = 'MCLE', pred_surv_args = list(),
                          mtryStart = NULL, stepFactor = 2,
                          time.eval = NULL, time.tau = NULL,
                          ntreeTry = 100L,
                          bootstrap = c("by.sub", "by.root", "none", "by.user"),
                          samptype = c("swor","swr"),
                          sampfrac = 0.632,
                          samp = NULL,
                          na.action = "na.omit",
                          trace = TRUE,
                          doBest = FALSE,
                          plot = FALSE,
                          applyfun = NULL, cores = NULL,
                          control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                            mincriterion = 0, saveinfo = FALSE,
                                                            minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                            minbucket = max(ceiling(sqrt(nrow(data))), 7))) {

  Call <- match.call()

  indx <- match(c('formula', 'id'), names(Call), nomatch = 0)
  if (indx[1] == 0) stop("a formula argument is required")

  yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]

  if (length(yvar.names) == 4){
    yvar.names = yvar.names[2:4]
  }
  n <- nrow(data)

  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)

  Ltimes <- data[, yvar.names[1]]
  Rtimes <- data[, yvar.names[2]]

  xvar.names <- attr(terms(formula), 'term.labels')
  nvar <- length(xvar.names)

  if (is.null(mtryStart)){
    mtryStart <- ceiling(sqrt(nvar))
  }

  if (indx[2] == 0){
    data$id <- 1:n
  } else {
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }

  data <- data[, c("id", yvar.names, xvar.names)]

  if (na.action == "na.omit") {
    takeid = which(complete.cases(data) == 1)
  } else if (na.action == "na.pass") {
    takeid = 1:n
  } else {
    stop("na.action can only be either 'na.omit' or 'na.pass'.")
  }

  id.sub <- unique(data$id[takeid])
  n.seu <- length(takeid)
  ## number of subjects
  n.sub <- length(id.sub)

  Rtimes <- Rtimes[takeid]
  if (n.seu == n.sub){
    if (is.null(time.eval)){
      time.eval <- c(0, sort(unique(Rtimes)))
    }
  } else {
    if (is.null(time.eval)){
      time.eval <- c(0, sort(unique(Rtimes)), seq(max(Rtimes), 1.5 * max(Rtimes), length.out = 50)[-1])
    }
    if (is.null(time.tau)){
      time.tau <- sapply(1:n.sub, function(ii){
        1.5 * max(Rtimes[data$id[takeid] == id.sub[ii]])
      })
    }
  }

  errorOOB_mtry <- function(eformula, edata, id,
                            eperm_test_est, eperm_test_args,
                            epred_surv_est, epred_surv_args,
                            emtryTest,
                            etpnt, etau,
                            entreeTry, econtrol,
                            ebootstrap,
                            esamptype,
                            esampfrac,
                            esamp,
                            ena.action, eapplyfun, ecores){
    cfOOB <- lbrccif(formula = eformula, data = edata, id = id,
                     perm_test_est = eperm_test_est, perm_test_args = eperm_test_args,
                     pred_surv_est = epred_surv_est, pred_surv_args = epred_surv_args,
                     mtry = emtryTest,
                     ntree = entreeTry,
                     control = econtrol,
                     bootstrap = ebootstrap,
                     samptype = esamptype,
                     sampfrac = esampfrac,
                     samp = esamp,
                     na.action = ena.action,
                     applyfun = eapplyfun,
                     cores = ecores)

    predOOB <- predictProb_LBRC(object = cfOOB, time.eval = etpnt, time.tau = etau, OOB = TRUE,
                                pred_surv_est = epred_surv_est, pred_surv_args = epred_surv_args)
    errorOOB <- sbrier_lbrc(obj = predOOB$survival.obj, id = predOOB$survival.id,
                            pred = predOOB, type = "IBS")
    rm(cfOOB)
    rm(predOOB)
    gc()
    return(errorOOB)
  }

  errorOld <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                            eperm_test_est = perm_test_est, eperm_test_args = perm_test_args,
                            epred_surv_est = pred_surv_est, epred_surv_args = pred_surv_args,
                            emtryTest = mtryStart,
                            etpnt = time.eval,
                            etau = time.tau,
                            entreeTry = ntreeTry,
                            econtrol = control,
                            ebootstrap = bootstrap,
                            esamptype = samptype,
                            esampfrac = sampfrac,
                            esamp = samp,
                            ena.action = na.action,
                            eapplyfun = applyfun,
                            ecores = cores)
  if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
  if (trace) {
    cat("mtry = ", mtryStart, " OOB Brier score = ",
        errorOld, "\n")
  }

  oobError <- list()
  oobError[[1]] <- errorOld
  names(oobError)[1] <- mtryStart

  for (direction in c("left", "right")) {
    if (trace) cat("Searching", direction, "...\n")
    mtryCur <- mtryStart
    while (mtryCur != nvar) {
      mtryOld <- mtryCur
      mtryCur <- if (direction == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(nvar, floor(mtryCur * stepFactor))
      }
      if (mtryCur == mtryOld) break

      errorCur <- errorOOB_mtry(eformula = formula, edata = data, id = id,
                                eperm_test_est = perm_test_est, eperm_test_args = perm_test_args,
                                epred_surv_est = pred_surv_est, epred_surv_args = pred_surv_args,
                                emtryTest = mtryCur,
                                etpnt = time.eval,
                                etau = time.tau,
                                entreeTry = ntreeTry,
                                econtrol = control,
                                ebootstrap = bootstrap,
                                esamptype = samptype,
                                esampfrac = sampfrac,
                                esamp = samp,
                                ena.action = na.action,
                                eapplyfun = applyfun,
                                ecores = cores)

      if (trace) {
        cat("mtry = ", mtryCur, "\tOOB error = ", errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      errorOld <- errorCur
    }
  }
  mtry <- sort(as.numeric(names(oobError)))
  res_all <- unlist(oobError[as.character(mtry)])
  res_all <- cbind(mtry = mtry, OOBError = res_all)
  res <- res_all[which.min(res_all[, 2]), 1]

  if (plot) {
    res <- res_all
    plot(res_all, xlab = expression(m[try]), ylab = "OOB Error", type = "o", log = "x", xaxt = "n")
    axis(1, at=res_all[, "mtry"])
  }

  if (doBest)
    res <- lbrccif(formula = formula, data = data, id = id,
                   perm_test_est = perm_test_est, perm_test_args = perm_test_args,
                   pred_surv_est = pred_surv_est, pred_surv_args = pred_surv_args,
                   mtry = res, ntree = ntreeTry,
                   control = control,
                   bootstrap = bootstrap,
                   samptype = samptype,
                   sampfrac = sampfrac,
                   samp = samp,
                   na.action = na.action,
                   applyfun = applyfun,
                   cores = cores)

  return(res)
}
