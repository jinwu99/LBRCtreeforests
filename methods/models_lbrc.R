library(partykit)
library(survival)

# logrank score for LTRC data
.logrank_trafo2 <- function(x2){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    # unique.times <- unique(x2[,2][which(x2[, 3] == 1)])
    # D <- rep(NA, length(unique.times))
    # R <- rep(NA, length(unique.times))
    # for(j in 1:length(unique.times)){
    #   D[j] = sum(unique.times[j] == x2[, 2])
    # }
    # for(k in 1:length(unique.times) ){
    #   value <- unique.times[k]
    #   R[k] <- sum(apply(x2[, 1:2], 1, function(interval){interval[1] < value & value <= interval[2]}))
    # }

    # faster with sapply
    delta <- x2[,3]
    A <- x2[,1]
    Z <- x2[,2]
    U <- unique(x2[,2][which(delta==1)])
    # number of failures at each distinct failure time
    D <- sapply(U, function(s) sum((Z==s)*delta))
    # number at risk at each distinct failure time
    R <- sapply(U, function(s) sum(A<=s & Z>=s))

    Ratio <- D / R

    Ratio <- Ratio[order(U)]
    Nelson.Aalen <- cumsum(Ratio)
    Event.time <- U[order(U)]
    Left <- sapply(x2[, 1], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
    Right <- sapply(x2[, 2], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})

    result<- x2[, 3] - (Right - Left)
  }
  return(as.double(result))
}

# score function of LBRC data with MCLE
.logrank_trafo3 <- function(x2){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    n <- dim(x2)[1]
    delta <- x2[,3]

    # unique 'failure' time
    U <- unique(x2[,2][which(delta==1)])
    A <- x2[,1]
    Z <- x2[,2]
    M <- max(Z)
    V <- Z-A

    # number of failures at each distinct failure time
    D <- sapply(U, function(s) sum((Z==s)*delta))
    # number at risk at each distinct failure time
    R <- sapply(U, function(s) sum((A<=s & Z>=s) + delta*(V<=s & Z>=s)))/2

    Ratio <- D / R
    Ratio <- Ratio[order(U)]
    Nelson.Aalen <- cumsum(Ratio)
    Event.time <- U[order(U)]
    logS_z <- -sapply(Z,function(t){if(t< min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
    result <- delta + logS_z
    # no need to compute (int S*logS / int S)
    # it cancels out in the standardization of test statistic
  }
  return(as.double(result))
}

# score function of LBRC data with MFLE
.logrank_trafo4 <- function(x2,eps=1e-8,max_iter=100){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    delta <- x2[,3]
    n <- length(delta)
    Z <- x2[,2]
    ix <- order(Z)

    res <- vardiCpp(x2,rep(1,n),eps=eps, max_iter=max_iter)
    t <- res$t
    p_order <- res$p
    lambda_order <- p_order / rev(cumsum(rev(p_order)))
    logS_z_order <- -cumsum(lambda_order)

    # return to original order
    logS_z <- logS_z_order[order(ix)]
    result <- delta + logS_z
    # no need to compute (int S*logS / int S)
    # it cancels out in the standardization of test statistic
  }
  return(as.double(result))
}

# Conditional inference trees/forests for LBRC (and LTRC) data
#
# Build a conditional inference tree (CIT) or forest (CIF) for length-biased
# right-censored (LBRC) data. For comparison/benchmarking, fitting on
# left-truncated right-censored (LTRC) data is also supported.
#
# Implementation adapted from LTRCforests:
# https://github.com/weichiyao/TimeVaryingData_LTRCforests/blob/main/pkg/LTRCforests/R/ltrccif.R
#
# Survival estimators:
# - "MCLE": nonparametric maximum composite conditional likelihood estimator of S(t).
# - "MFLE": nonparametric maximum full-likelihood estimator of the unbiased survival function.
# - "KM":   risk-adjusted Kaplan–Meier (default in LTRC-CIT/CIF).
#
# Key arguments:
# - pred_surv_est: Character. Estimator for unbiased survival prediction;
#   one of c("MCLE", "MFLE", "KM").
# - pred_surv_args: List. EM routine controls when pred_surv_est == "MFLE"
#   (ignored otherwise). Fields: eps (convergence tolerance),
#   max_iter (maximum EM iterations).
#
# Details:
# - The algorithm uses permutation tests at internal nodes with the chosen survival estimator
#   and predicts unbiased survival function using the estimator specified for prediction.
# - For other arguments, see LTRCforests::ltrccif.
lbrccit <- function(formula, data, id,
                    perm_test_est = "MCLE", perm_test_args = list(),
                    na.action = "na.omit",
                    control = partykit::ctree_control()){
  Call <- match.call()
  Call[[1]] <- as.name('lbrccit')
  indx <- match(c('formula', 'data', 'id'),
                names(Call), nomatch = 0L)
  if (indx[1] == 0) stop("a formula argument is required")
  Call$formula <- eval(formula)

  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)

  mf <- eval.parent(temp)
  y <- model.extract(mf, 'response')
  if (!is.Surv(y)) stop("Response must be a survival object")
  if (!attr(y, "type") == "counting") stop("The Surv object must be of type 'counting'.")
  rm(y)

  # pull y-variable names
  yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]
  if (length(yvar.names) == 4){
    yvar.names = yvar.names[2:4]
  }

  Status <- data[, yvar.names[3]]
  if (sum(Status) == 0) stop("All observations are right-censored with event = 0!")

  n <- nrow(data)
  if (indx[3] == 0){
    ## If id is not present, then we add one more variable
    data$id <- 1:n
  } else {
    ## If id is present, then we rename the column to be id
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }

  # extract the x-variable names
  xvar.names <- attr(terms(formula), 'term.labels')
  rm(temp)
  data <- data[, c("id", yvar.names, xvar.names)]

  # permutation test statistics for unbiased variable selection
  if(perm_test_est == 'MFLE'){
    logrank_trafo <- function(...) do.call(.logrank_trafo4, c(list(...),perm_test_args))
  }else if(perm_test_est == "KM"){
    logrank_trafo <- .logrank_trafo2
  }else{ # MCLE
    logrank_trafo <- .logrank_trafo3
  }

  h <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, eps = 1e-7, max_iter = 100, ...) {
    if (is.null(weights)) weights <- rep(1, NROW(y))
    s <- logrank_trafo(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  ret <- partykit::ctree(formula = formula, data = data, na.action = na.action,
                         control = control, ytrafo = h)

  ret$formulaLBRC <- formula
  ret$info$call <- Call
  ret$perm_test_est <- perm_test_est
  if (na.action == "na.omit"){
    ret$data$id <- data$id[complete.cases(data) == 1]
  } else {
    ret$data$id <- data$id
  }
  class(ret) <- c("lbrccit",class(ret))
  ret
}

lbrccif <- function(formula, data, id,
                    perm_test_est = "MCLE", perm_test_args = list(),
                    pred_surv_est = "MCLE", pred_surv_args = list(),
                    mtry = NULL, ntree = 100L,
                    bootstrap = c("by.sub","by.root","by.user","none"),
                    samptype = c("swor","swr"),
                    sampfrac = 0.632,
                    samp = NULL,
                    na.action = "na.omit",
                    stepFactor = 2,
                    trace = TRUE,
                    applyfun = NULL, cores = NULL,
                    control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                                      minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                      minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                                      mincriterion = 0, saveinfo = FALSE)){
  Call <- match.call()
  Call[[1]] <- as.name('lbrccif')
  indx <- match(c('formula', 'data', 'id'),
                names(Call), nomatch = 0L)
  if (indx[1] == 0) stop("a formula argument is required")
  Call$formula <- eval(formula)

  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)

  mf <- eval.parent(temp)
  y <- model.extract(mf, 'response')
  if (!is.Surv(y)) stop("Response must be a survival object")
  if (!attr(y, "type") == "counting") stop("The Surv object must be of type 'counting'.")
  rm(y)

  yvar.names <- all.vars(formula(paste(as.character(formula)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]

  if (length(yvar.names) == 4){
    yvar.names = yvar.names[2:4]
  }

  Status <- data[, yvar.names[3]]
  if (sum(Status) == 0) stop("All observations are right-censored with event = 0!")

  n <- nrow(data)

  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)

  if (indx[3] == 0){
    data$id <- 1:n
  } else {
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }

  xvar.names <- attr(terms(formula), 'term.labels')
  rm(temp)
  data <- data[, c("id", yvar.names, xvar.names)]

  if (length(data$id) == length(unique(data$id))){
    if (bootstrap == "by.sub") bootstrap <- "by.root"
  } else {
    id.sub <- unique(data$id)
    n.sub <- length(id.sub)
  }

  if (samptype == "swor"){
    perturb = list(replace = FALSE, fraction = sampfrac)
  } else if (samptype == "swr"){
    perturb = list(replace = TRUE)
  } else {
    stop("samptype must set to be either 'swor' or 'swr'\n")
  }

  if (bootstrap == "by.sub"){
    size <- n.sub
    if (!perturb$replace) size <- floor(n.sub * perturb$fraction)
    samp <- replicate(ntree,
                      sample(id.sub, size = size,
                             replace = perturb$replace),
                      simplify = FALSE) # a list of length ntree
    samp <- lapply(samp, function(y) unlist(sapply(y, function(x) which(data$id %in% x), simplify = FALSE)))
    samp <- sapply(samp, function(x) as.integer(tabulate(x, nbins = n))) # n x ntree
  } else if (bootstrap == "none"){
    samp <- matrix(1, nrow = n, ncol = ntree)
  } else if (bootstrap == "by.user") {
    if (is.null(samp)) {
      stop("samp must not be NULL when bootstrapping by user\n")
    }
    if (is.matrix(samp)){
      if (!is.matrix(samp)) stop("samp must be a matrx\n")
      if (any(!is.finite(samp))) stop("samp must be finite\n")
      if (any(samp < 0)) stop("samp must be non-negative\n")
      if (all(dim(samp) != c(n, ntree))) stop("dimension of samp must be n x ntree\n")
      samp <- as.matrix(samp)  # transform into matrix
    }
  } else if (bootstrap == "by.root"){
    samp <- rep(1, n)
  } else {
    stop("Wrong bootstrap is given!\n ")
  }

  if (is.null(mtry)){
    mtry <- tune.lbrccif(formula = formula, data = data, id = id,
                         perm_test_est = perm_test_est, perm_test_args = perm_test_args,
                         pred_surv_est = pred_surv_est, pred_surv_args = pred_surv_args,
                         control = control, ntreeTry = ntree,
                         bootstrap = "by.user",
                         samptype = samptype,
                         sampfrac = sampfrac,
                         samp = samp,
                         na.action = na.action,
                         stepFactor = stepFactor,
                         applyfun = applyfun,
                         cores = cores,
                         trace = trace)
    print(sprintf("mtry is tuned to be %1.0f", mtry))
  }

  # permutation test statistics for unbiased variable selection
  if(perm_test_est == 'MFLE'){
    logrank_trafo <- function(...) do.call(.logrank_trafo4, c(list(...),perm_test_args))
  }else if(perm_test_est == "KM"){
    logrank_trafo <- .logrank_trafo2
  }else{ # MCLE
    logrank_trafo <- .logrank_trafo3
  }

  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (all(is.na(weights)) == 1) weights <- rep(1, NROW(y))
    s <- logrank_trafo(y[weights > 0, , drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  ret <- partykit::cforest(formula, data,
                           weights = samp,
                           perturb = perturb,
                           ytrafo = h2,
                           control = control,
                           na.action = na.action,
                           mtry = mtry,
                           ntree = ntree,
                           applyfun = applyfun,
                           cores = cores)
  ret$mtry <- mtry
  ret$formulaLBRC <- formula
  ret$info$call <- Call
  ret$info$bootstrap <- bootstrap
  ret$info$samptype <- samptype
  ret$info$sampfrac <- sampfrac
  ret$perm_test_est <- perm_test_est
  ret$tune_pred_surv_est <- ifelse(is.null(mtry), NULL, pred_surv_est)
  if (na.action == "na.omit"){
    ret$data$id <- data$id[complete.cases(data) == 1]
  } else {
    ret$data$id <- data$id
  }
  class(ret) <- "lbrccif"
  ret
}



