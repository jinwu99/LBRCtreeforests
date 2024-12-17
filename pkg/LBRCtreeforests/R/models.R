.logrank_trafo2 <- function(x2){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    unique.times <- unique(x2[,2][which(x2[, 3] == 1)])

    # D <- rep(NA, length(unique.times))
    # R <- rep(NA, length(unique.times))
    #
    # for(j in 1:length(unique.times)){
    #   D[j] = sum(unique.times[j] == x2[, 2])
    # }
    #
    # for(k in 1:length(unique.times) ){
    #   value <- unique.times[k]
    #   R[k] <- sum(apply(x2[, 1:2], 1, function(interval){interval[1] < value & value <= interval[2]}))
    # }

    delta <- x2[,3]
    A <- x2[,1]
    Z <- x2[,2]
    U <- unique(x2[,2][which(delta==1)])
    # number of failures at each distinct failure time
    D <- sapply(U, function(s) sum((Z==s)*delta))
    # number at risk at each distinct failure time
    R <- sapply(U, function(s) sum(A<=s & Z>=s))

    Ratio <- D / R

    Ratio <- Ratio[order(unique.times)]
    Nelson.Aalen <- cumsum(Ratio)
    Event.time <- unique.times[order(unique.times)]
    Left <- sapply(x2[, 1], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
    Right <- sapply(x2[, 2], function(t){if(t < min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})

    result<- x2[, 3] - (Right - Left)
  }
  return(as.double(result))
}

.logrank_trafo3 <- function(x2){
  if (sum(x2[, 3] == 1) == 0) {
    result <- x2[,3]
  } else {
    n <- dim(x2)[1]
    delta <- x2[,3]

    U <- unique(x2[,2][which(delta==1)]) # unique 'failure' time
    # M <- max(U)
    A <- x2[,1]
    Z <- x2[,2]
    M <- max(Z)
    V <- Z-A
    AV <- sort(unique(c(0,A,V)))

    # number of failures at each distinct failure time
    dN <- sapply(U, function(s) sum((Z==s)*delta))
    # number at risk at each distinct failure time
    R <- sapply(U, function(s) sum(A<=s & Z>=s))

    dQ <- sapply(AV, function(s) sum(A==s) + sum((V==s)*delta))
    K <- sapply(AV, function(s) sum(A>=s) + sum(V>=s))
    K[K<=0] <- eps
    S_A <- sapply(U, function(s){
      prod(1-dQ[AV<=s]/K[AV<=s])
    })
    R2 <- sapply(U, function(s) sum(Z>=s)) - S_A*n
    R2[R2<=0] <- eps

    Ratio <- dN/R2
    Nelson.Aalen <- cumsum(Ratio)
    Event.time <- U[order(U)]

    logS_z <- -sapply(Z,function(t){if(t< min(Event.time)) return(0) else return(Nelson.Aalen[max(which(Event.time <= t))])})
    result <- delta + logS_z
  }
  return(as.double(result))
}

#' Fit a LBRC conditional inference tree and forest
#' Fitting LTRC conditional inference tree and forest is also available
#' We have made used of codes in LTRCforests package, adding some additional functions appropriate for LBRC data
#' @import partykit
#' @import survival
#' @import LTRCforests
#' @export
ltrccit <- function(formula, data, id, lenbias = FALSE, na.action = "na.omit",
                    control = partykit::ctree_control()){
  #requireNamespace("inum")

  Call <- match.call()
  Call[[1]] <- as.name('ltrccit')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
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
  # Times <- data[, yvar.names[2]]
  if (sum(Status) == 0) stop("All observations are right-censored with event = 0!")

  n <- nrow(data)

  if (indx[3] == 0){
    ## If id is not present, then we add one more variable
    # mf$id <- 1:nrow(mf) ## Relabel
    data$id <- 1:n
  } else {
    ## If id is present, then we rename the column to be id
    # names(mf)[names(mf) == "(id)"] <- "id"
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }

  # extract the x-variable names
  xvar.names <- attr(terms(formula), 'term.labels')
  rm(temp)
  data <- data[, c("id", yvar.names, xvar.names)]

  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (is.null(weights)) weights <- rep(1, NROW(y))
    s <- .logrank_trafo2(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  h3 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (is.null(weights)) weights <- rep(1, NROW(y))
    s <- .logrank_trafo3(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  if(lenbias == TRUE){
    ret <- partykit::ctree(formula = formula, data = data,
                           na.action = na.action, ytrafo = h3, control = control)
  }else{
    ret <- partykit::ctree(formula = formula, data = data,
                           na.action = na.action, ytrafo = h2, control = control)
  }

  ret$formulaLTRC <- formula
  ret$info$call <- Call
  ret$lenbias <- lenbias
  if (na.action == "na.omit"){
    ret$data$id <- data$id[complete.cases(data) == 1]
  } else {
    ret$data$id <- data$id
  }
  class(ret) <- c("ltrccit",class(ret))
  ret
}

ltrccif <- function(formula, data, id, lenbias = FALSE,
                    FUN = .pred_Surv_nolog,
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
                                                      # minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                                      # minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                                      minsplit = nrow(data) * 0.15,
                                                      minbucket = nrow(data) * 0.06,
                                                      minprob = 0.01,
                                                      # maxdepth = 2,
                                                      # mincriterion = 0,
                                                      saveinfo = FALSE)){

  #requireNamespace("inum")

  Call <- match.call()
  Call[[1]] <- as.name('ltrccif')  #make nicer printout for the user
  # create a copy of the call that has only the arguments we want,
  #  and use it to call model.frame()
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
  # Times <- data[, yvar.names[2]]
  if (sum(Status) == 0) stop("All observations are right-censored with event = 0!")

  n <- nrow(data)

  ## if not specified, the first one will be used as default
  bootstrap <- match.arg(bootstrap)
  samptype <- match.arg(samptype)

  ## The following code to define id does not work since it could not handle missing values
  # id <- model.extract(mf, 'id')

  # this is a must, otherwise id cannot be passed to the next level in tune.ltrccif
  if (indx[3] == 0){
    ## If id is not present, then we add one more variable
    # mf$id <- 1:nrow(mf) ## Relabel
    data$id <- 1:n
  } else {
    ## If id is present, then we rename the column to be id
    # names(mf)[names(mf) == "(id)"] <- "id"
    names(data)[names(data) == deparse(substitute(id))] <- "id"
  }

  # extract the x-variable names
  xvar.names <- attr(terms(formula), 'term.labels')
  rm(temp)
  data <- data[, c("id", yvar.names, xvar.names)]

  ## bootstrap case
  if (length(data$id) == length(unique(data$id))){ # time-invariant LTRC data
    # it includes the case 1) when id = NULL, which is that id is not specified
    #                      2) when id is specified, but indeed LTRC time-invariant
    if (bootstrap == "by.sub") bootstrap <- "by.root"
  } else { # time-varying subject data
    if (lenbias == TRUE) stop("Time invariant data is only allowed for LB model")
    id.sub <- unique(data$id)
    ## number of subjects
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
    mtry <- tune.ltrccif(formula = formula, data = data, id = id,
                         lenbias = lenbias,
                         FUN = FUN,
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

  # for left truncated data
  h2 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (all(is.na(weights)) == 1) weights <- rep(1, NROW(y))
    s <- .logrank_trafo2(y[weights > 0, , drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  # for length biased data
  h3 <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (all(is.na(weights)) == 1) weights <- rep(1, NROW(y))
    s <- .logrank_trafo3(y[weights > 0, , drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  if (lenbias == TRUE){
    ret <- partykit::cforest(formula, data,
                             weights = samp,
                             perturb = perturb,
                             ytrafo = h3,
                             control = control,
                             na.action = na.action,
                             mtry = mtry,
                             ntree = ntree,
                             applyfun = applyfun,
                             cores = cores)
  }else{
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
  }

  ret$formulaLTRC <- formula
  ret$info$call <- Call
  ret$info$bootstrap <- bootstrap
  ret$info$samptype <- samptype
  ret$info$sampfrac <- sampfrac
  ret$lenbias <- lenbias
  if (na.action == "na.omit"){
    ret$data$id <- data$id[complete.cases(data) == 1]
  } else {
    ret$data$id <- data$id
  }
  class(ret) <- "ltrccif"
  ret
}

