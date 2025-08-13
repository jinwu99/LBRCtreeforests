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
#
# For explanations of other arguments, see LTRCforests::ltrccif.
predictProb_LBRC <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                             pred_surv_est = "MCLE", pred_surv_args = list(),
                             time.eval, time.tau = NULL){
  UseMethod("predictProb_LBRC", object)
}

predictProb_LBRC.lbrccit <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                     pred_surv_est = "MCLE", pred_surv_args = list(),
                                     time.eval, time.tau = NULL){

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

  rm(object)


  N <- length(unique(newdata[, "id"]))

  if (is.null(time.tau)){
    time.tau <- rep(max(time.eval), N)
  } else {
    if (N != length(time.tau)) stop("time.tau should be a vector of length equaling to number of SUBJECT observation!")
  }

  Shat <- sapply(1:N, function(Ni) .shatfunc(Ni, data = newdata, pred = pred, tpnt = time.eval, tau = time.tau))
  obj <- Surv(newdata[, yvar.names[1]],
              newdata[, yvar.names[2]],
              newdata[, yvar.names[3]])
  RES <- list(survival.probs = Shat,
              survival.times = time.eval,
              survival.tau = time.tau,
              survival.obj = obj,
              survival.id = newdata$id,
              survival.est = pred_surv_est)
  rm(newdata)
  rm(Shat)
  rm(time.eval)
  rm(time.tau)
  return(RES)
}

predictProb_LBRC.lbrccif <- function(object, newdata = NULL, newdata.id, OOB = FALSE,
                                     pred_surv_est = "MCLE", pred_surv_args = list(),
                                     time.eval, time.tau = NULL){
  if(pred_surv_est == "MFLE"){
    pred_FUN <- function(...) do.call(.pred_Surv_LBRC_MFLE, c(list(...), pred_surv_args))
  }else if(pred_surv_est == "KM"){
    pred_FUN <- .pred_Surv_nolog
  }else if(pred_surv_est == "MCLE"){
    pred_FUN <- .pred_Surv_LBRC_MCLE
  }

  if(OOB){
    m <- as.matrix(object$data[[1]])
    # Rebuild the Surv(start, stop, event) with start = 0
    # Needed to compute the OOB unbiased survival function during tuning
    object$data[[1]] <- Surv(time  = rep(0, nrow(m)),
                             time2 = m[,2],
                             event = m[,3])
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

  Shat <- sapply(1:N, function(Ni) .shatfunc(Ni, data = newdata, pred = pred, tpnt = time.eval, tau = time.tau))
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
  if(nj == 1){
    II <- which(id.seu == id.sub[Ni])[jall[nj]]
    Shat_i = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[nj]])
    Shat_i[1] <- ifelse(Shat_i[1]<=0, 1e-10, Shat_i[1])
    Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
  } else if (nj > 1) {
    ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
    ShatR_temp[1, 1] <- 1
    qL = rep(0, nj)
    for (j in 1:nj){
      II <- which(id.seu == id.sub[Ni])[1] + jall[j] - 1
      Shat_j = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[j]])
      qL[j] <- Shat_j[1]
      jR = ipred::getsurv(pred[[II]], TestT[j + 1])
      ShatR_temp[1, j + 1] = jR / qL[j]
      Shat_temp[1, r.ID == jall[j]] <- Shat_j / qL[j]
    }

    ql0 <- which(qL == 0)
    if (length(ql0) > 0){
      if (any(qL > 0)){
        maxqlnot0 <- max(which(qL > 0))

        ql0lmax <- ql0[ql0 < maxqlnot0]
        ql0mmax <- ql0[ql0 >= maxqlnot0]
        ShatR_temp[1, ql0lmax + 1] <- 1
        Shat_temp[1, r.ID %in% jall[ql0lmax]] <- 1
        ShatR_temp[1, ql0mmax + 1] <- 0
        Shat_temp[1, r.ID %in% jall[ql0mmax]] <- 0
      } else {
        ShatR_b[1, 2:(nj + 1)] <- 0
        Shat_temp[1, r.ID %in% jall] <- 0
      }
    }
    m <- cumprod(ShatR_temp[1, 1:nj])
    for (j in 1:nj){
      Shat_temp[1, r.ID == jall[j]] <- Shat_temp[1, r.ID == jall[j]] * m[j]
    }
  }

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


.pred_Surv_LBRC_MFLE <- function(y,w,eps=1e-7,max_iter=100){
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
