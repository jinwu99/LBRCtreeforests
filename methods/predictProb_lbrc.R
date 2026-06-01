eps <- 1e-10

# Compute a survival curve from lbrrcit and lbrccif models
#
# Implementation adapted from LTRCforests:
# https://github.com/weichiyao/TimeVaryingData_LTRCforests/blob/main/pkg/LTRCforests/R/predictProb.R
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
# - target: Character. Survival function target for prediction;
#   one of c("observed", "unbiased"). "observed" returns the conditional
#   survival function given delayed entry, P(T > t | T > A), whereas
#   "unbiased" returns the underlying survival function S(t).
#
# For explanations of other arguments, see LTRCforests::ltrccif.
predictProb_LBRC <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                             pred_surv_est = "MCLE", pred_surv_args = list(),
                             time.eval, time.tau = NULL,
                             target = c("observed", "unbiased")){
  UseMethod("predictProb_LBRC", object)
}

predictProb_LBRC.lbrccit <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                     pred_surv_est = "MCLE", pred_surv_args = list(),
                                     time.eval, time.tau = NULL,
                                     target = c("observed", "unbiased")){

  if(pred_surv_est == "MFLE"){
    pred_FUN <- function(...) do.call(.pred_Surv_LBRC_MFLE, c(list(...), pred_surv_args))
  }else if(pred_surv_est == "KM"){
    pred_FUN <- .pred_Surv_nolog
  }else if(pred_surv_est == "MCLE"){
    pred_FUN <- .pred_Surv_LBRC_MCLE
  }

  pred <- predict(object, type='prob', newdata=newdata, FUN = pred_FUN)

  xvar.names <- attr(object$terms,"term.labels")
  yvar.names <- as.character(object$formulaLBRC[[2]])[2:4]
  idname <- "id"

  if (is.null(newdata)){
    newdata <- as.data.frame(as.matrix(object$data[, c(1, ncol(object$data)), drop = FALSE]))
    names(newdata) = c(yvar.names, idname)
  } else {
    if (missing(newdata.id)){
      newdata$id <- 1:nrow(newdata)
    } else {
      names(newdata)[names(newdata) == deparse(substitute(newdata.id))] <- idname
    }
    newdata <- as.data.frame(newdata[, c(yvar.names, idname)])
  }

  target <- match.arg(target)
  preddata <- newdata
  if(target == "unbiased"){
    # Target unbiased survival function by forcing all truncation time to 0
    preddata[, yvar.names[1]] <- 0
  }

  rm(object)

  N <- length(unique(newdata[, "id"]))

  if (is.null(time.tau)){
    time.tau <- rep(max(time.eval), N)
  } else {
    if (N != length(time.tau)) stop("time.tau should be a vector of length equaling to number of SUBJECT observation!")
  }

  Shat <- sapply(1:N, function(Ni) .shatfunc(Ni, data = preddata, pred = pred, tpnt = time.eval, tau = time.tau))
  obj <- Surv(newdata[, yvar.names[1]],
              newdata[, yvar.names[2]],
              newdata[, yvar.names[3]])
  RES <- list(survival.probs = Shat,
              survival.times = time.eval,
              survival.tau = time.tau,
              survival.obj = obj,
              survival.id = newdata$id,
              survival.est = pred_surv_est,
              target = target)
  rm(newdata)
  rm(Shat)
  rm(time.eval)
  rm(time.tau)
  return(RES)
}

predictProb_LBRC.lbrccif <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                     pred_surv_est = "MCLE", pred_surv_args = list(),
                                     time.eval, time.tau = NULL,
                                     target = c("observed", "unbiased")){
  if(pred_surv_est == "MFLE"){
    pred_FUN <- function(...) do.call(.pred_Surv_LBRC_MFLE, c(list(...), pred_surv_args))
  }else if(pred_surv_est == "KM"){
    pred_FUN <- .pred_Surv_nolog
  }else if(pred_surv_est == "MCLE"){
    pred_FUN <- .pred_Surv_LBRC_MCLE
  }

  pred <- partykit::predict.cforest(object = object, newdata = newdata, OOB = OOB, type = "prob",
                                    FUN = pred_FUN)
  xvar.names <- attr(object$terms,"term.labels")
  yvar.names <- as.character(object$formulaLBRC[[2]])[2:4]
  idname <- "id"

  if (is.null(newdata) || OOB){
    newdata <- as.data.frame(as.matrix(object$data[, c(1, ncol(object$data)), drop = FALSE]))
    names(newdata) = c(yvar.names, idname)
  } else {
    if (missing(newdata.id)){
      newdata$id <- 1:nrow(newdata)
    } else {
      names(newdata)[names(newdata) == deparse(substitute(newdata.id))] <- idname
    }
    newdata <- as.data.frame(newdata[, c(yvar.names, idname)])
  }
  rm(object)

  N <- length(unique(newdata[, "id"]))

  if (is.null(time.tau)){
    time.tau <- rep(max(time.eval), N)
  } else {
    if (N != length(time.tau)) stop("time.tau should be a vector of length equaling to number of SUBJECT observation! \n
                                     In the time-varying case, check whether newdata.id has been correctly specified!")
  }

  target <- match.arg(target)
  preddata <- newdata
  if(target == "unbiased"){
    # Target unbiased survival function by forcing all truncation time to 0
    preddata[, yvar.names[1]] <- 0
  }

  Shat <- sapply(1:N, function(Ni) .shatfunc(Ni, data = preddata, pred = pred, tpnt = time.eval, tau = time.tau))
  obj <- Surv(newdata[, yvar.names[1]],
              newdata[, yvar.names[2]],
              newdata[, yvar.names[3]])
  RES <- list(survival.probs = Shat,
              survival.times = time.eval,
              survival.tau = time.tau,
              survival.obj = obj,
              survival.id = newdata$id)
  rm(newdata)
  rm(Shat)
  rm(time.eval)
  rm(time.tau)
  return(RES)
}

batched_predictProb_LBRC <- function(object,
                                     newdata,
                                     time.eval,
                                     time.tau,
                                     pred_surv_est,
                                     pred_surv_args = list(),
                                     batch_size = 50) {
  n_test <- nrow(newdata)
  num_batches <- ceiling(n_test / batch_size)
  all_survival_probs <- vector("list", num_batches)

  for (b in 1:num_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, n_test)
    chunk <- newdata[start_idx:end_idx, , drop = FALSE]

    pred_chunk <- predictProb_LBRC(object = object,
                                   newdata = chunk,
                                   time.eval = time.eval,
                                   time.tau = time.tau,
                                   pred_surv_est = pred_surv_est,
                                   pred_surv_args = pred_surv_args,
                                   target = "observed")
    all_survival_probs[[b]] <- pred_chunk$survival.probs
    rm(pred_chunk, chunk)
    gc()
  }

  combined_probs <- do.call(cbind, all_survival_probs)
  return(list(survival.probs = combined_probs))
}

# Edited LTRCforests::.shatfunc to compute Shat_i safely (prevents division by zero)
.shatfunc <- function(Ni, data, pred, tpnt, tau){
  id.seu <- data[, "id"]
  id.sub <- unique(id.seu)

  TestData <- data[id.seu == id.sub[Ni], ]

  TestT <- c(TestData[1, 1], TestData[, 2])
  TestTIntN <- nrow(TestData)

  tpnt <- tpnt[tpnt <= tau[Ni]]

  tpntL <- c(TestT, tpnt)
  torder <- order(tpntL)
  tpntLod <- tpntL[torder]
  tlen <- length(tpntLod)

  Shat_temp <- matrix(0, nrow = 1, ncol = tlen)

  r.ID <- findInterval(tpntLod, TestT)
  r.ID[r.ID > TestTIntN] <- TestTIntN

  jall <- unique(r.ID[r.ID > 0])
  nj <- length(jall)

  Shat_temp[1, r.ID == 0] <- 1
  II <- which(id.seu == id.sub[Ni])[jall[nj]]
  Shat_i = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[nj]])
  Shat_i[1] <- ifelse(Shat_i[1]<=0, 1e-10, Shat_i[1])
  Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]

  return(Shat_temp[1, -match(TestT, tpntLod)])
  rm(Shat_temp)
  rm(ShatR_temp)
  rm(id.seu)
  rm(id.sub)
}


.pred_Surv_nolog <- function(y, w) {
  if (length(y) == 0) return(NA)
  survfit(y ~ 1, weights = w, subset = w > 0, conf.type = "none", se.fit = FALSE)
}


.pred_Surv_LBRC_MCLE <- function(y,w){
  if (length(y) == 0) return(NA)
  idx = which(w>0)
  y = y[idx,]
  w = w[idx]

  n <- sum(w)
  delta <- y[,3]

  if(sum(delta)==0){
    S_pred <- as.double(rep(1,dim(y)[1]))
    return(structure(
      list(time=y[,2], surv=S_pred),
      class = c("survfit_lb", "survfit")
    ))
  }

  U <- unique(y[,2][which(delta==1)])
  A <- y[,1]
  Z <- y[,2]
  V <- Z-A

  dN <- sapply(U, function(s) sum(w*(Z==s)*delta))
  R <- sapply(U, function(s) sum( w*( (A<=s & Z>=s) + delta*(V<=s & Z>=s) ) ) )/2

  dN_R <- dN/R
  dN_R[dN_R>1] <- 1

  # with unique event times
  Y <- sort(unique(Z))
  S_pred <- sapply(Y, function(x){
    prod(1-dN_R[U<=x])
  })

  S_pred[S_pred<0] <- 0; S_pred[S_pred>1] <- 1

  return(structure(
    list(time=Y, surv=S_pred),
    class = c("survfit_lb", "survfit")
  ))
}


.pred_Surv_LBRC_MFLE <- function(y,w,eps=1e-7,max_iter=30){
  if (length(y) == 0) return(NA)
  idx = which(w>0)
  y = y[idx,]
  w = w[idx]

  n <- sum(w)
  delta <- y[,3]

  if(sum(delta)==0){
    S_pred <- as.double(rep(1,dim(y)[1]))
    return(structure(
      list(time=y[,2], surv=S_pred),
      class = c("survfit_lb", "survfit")
    ))
  }

  Z <- y[,2]
  Z_sort <- sort(Z)

  res <- vardiCpp(y,w,eps = eps,max_iter = max_iter)
  Y <- res$t
  S_pred <- res$S
  S_pred[S_pred<0] <- 0; S_pred[S_pred>1] <- 1

  return(structure(
    list(time=Y, surv=S_pred),
    class = c("survfit_lb", "survfit")
  ))
}


print.survfit_lb <- function(x, ...){
  cat("LB Survival Object\n")
  cat("Median survival time :",format(median(x$time),digits=3))
}
