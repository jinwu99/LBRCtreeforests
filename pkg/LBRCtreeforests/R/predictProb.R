eps <- 1e-10

# modification of ipred::getsurv
getsurv <- function (obj, times){
  if (!inherits(obj, c("survfit", "survfit_lb")))
    stop("obj is not of class survfit or LBsurv")
  class(obj) <- NULL
  lt <- length(times)
  nsurv <- times
  if (length(times) == length(obj$time)) {
    if (all(times == obj$time))
      return(obj$surv)
  }
  inside <- times %in% obj$time
  for (i in (1:lt)) {
    if (inside[i])
      nsurv[i] <- obj$surv[obj$time == times[i]]
    else {
      less <- obj$time[obj$time < times[i]]
      if (length(less) == 0)
        nsurv[i] <- 1
      else nsurv[i] <- obj$surv[obj$time == max(less)]
    }
  }
  nsurv
}

#' Compute a Survival Curve from ltrrcit and ltrccif model
#' @import partykit
#' @import survival
#' @export
predictProb <- function(object, newdata = NULL, newdata.id,
                        FUN = .pred_Surv_nolog,
                        OOB = FALSE, time.eval, time.tau = NULL){
  UseMethod("predictProb", object)
}
#' @export
predictProb.ltrccit <- function(object, newdata = NULL, newdata.id,
                                FUN = .pred_Surv_nolog,
                                OOB = FALSE, time.eval, time.tau = NULL){
  pred <- predict(object, type='prob', newdata=newdata, FUN=FUN)
  xvar.names <- attr(object$terms,"term.labels")
  yvar.names <- as.character(object$formulaLTRC[[2]])[2:4]
  idname <- "id"

  # missing values can be present in the prediction
  if (is.null(newdata)){
    # first column: Surv(tleft,tright,event), second column: (id)
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


  N <- length(unique(newdata[, "id"])) # number of subjects

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
              survival.id = newdata$id)
  rm(newdata)
  rm(Shat)
  rm(time.eval)
  rm(time.tau)
  return(RES)
}
#' @export
predictProb.ltrccif <- function(object, newdata = NULL, newdata.id,
                                FUN = .pred_Surv_nolog,
                                OOB = FALSE, time.eval, time.tau = NULL){
  pred <- partykit::predict.cforest(object = object, newdata = newdata, OOB = OOB, type = "prob",
                                    FUN = FUN)
  xvar.names <- attr(object$terms,"term.labels")
  yvar.names <- as.character(object$formulaLTRC[[2]])[2:4]
  idname <- "id"

  # missing values can be present in the prediction
  if (is.null(newdata) || OOB){
    # first column: Surv(tleft,tright,event), second column: (id)
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


  N <- length(unique(newdata[, "id"])) # number of subjects

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
# Edited on obtaining Shat_i to prevent of dividing 0 value
# this frequently happens for LBRC data
.shatfunc <- function(Ni, data, pred, tpnt, tau){
  ## This function is to compute the estimated survival probability of the Ni-th subject
  id.seu <- data[, "id"] # id
  id.sub <- unique(id.seu)

  ## the i-th data
  TestData <- data[id.seu == id.sub[Ni], ]

  TestT <- c(TestData[1, 1], TestData[, 2])
  TestTIntN <- nrow(TestData)

  tpnt <- tpnt[tpnt <= tau[Ni]]

  ################ Changes at July 29th
  tpntL <- c(TestT, tpnt)
  torder <- order(tpntL)
  tpntLod <- tpntL[torder]
  tlen <- length(tpntLod)

  ## Compute the estimated survival probability of the Ni-th subject
  Shat_temp <- matrix(0, nrow = 1, ncol = tlen)

  r.ID <- findInterval(tpntLod, TestT)
  r.ID[r.ID > TestTIntN] <- TestTIntN

  jall <- unique(r.ID[r.ID > 0])
  nj <- length(jall)

  ## Deal with left-truncation
  Shat_temp[1, r.ID == 0] <- 1
  if(nj == 1){
    ## Get the index of the Pred to compute Shat
    II <- which(id.seu == id.sub[Ni])[jall[nj]]
    Shat_i = getsurv(pred[[II]], tpntLod[r.ID == jall[nj]])
    # Survival function at truncation time (Shat_i[1]) can be 0
    Shat_i[1] <- ifelse(Shat_i[1]<=0, 1e-10, Shat_i[1])
    Shat_temp[1, r.ID == jall[nj]] <- Shat_i / Shat_i[1]
  } else if (nj > 1) {
    # c(1, S_{1}(R_{1}), ..., S_{nj}(R_{nj}))
    ShatR_temp <- matrix(0, nrow = 1, ncol = nj + 1)
    ShatR_temp[1, 1] <- 1

    # S_1(L_1), S_2(L_2), S_3(L_3), ..., S_{nj}(L_{nj})
    qL = rep(0, nj)
    for (j in 1:nj){
      ## Get the index of the Pred to compute Shat
      II <- which(id.seu == id.sub[Ni])[1] + jall[j] - 1
      Shat_j = ipred::getsurv(pred[[II]], tpntLod[r.ID == jall[j]])

      qL[j] <- Shat_j[1]
      # S_{j}(R_{j}), j=1,...nj-1
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
        # for(j in ql0){
        #   if (j < maxqlnot0) {
        #     ShatR_temp[1, j + 1] <- 1
        #     Shat_temp[1, r.ID == jall[j]] <- 1
        #   } else{
        #     ShatR_temp[1, j + 1] <- 0
        #     Shat_temp[1, r.ID == jall[j]] <- 0
        #   }
        # }
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

  # since: tpntLod[torder == 1] == TestData[1, 1]
  return(Shat_temp[1, -match(TestT, tpntLod)])
  rm(Shat_temp)
  rm(ShatR_temp)
  rm(id.seu)
  rm(id.sub)
}
#' @export
.pred_Surv_nolog <- function(y, w) {
  if (length(y) == 0) return(NA)
  survfit(y ~ 1, weights = w, subset = w > 0, conf.type = "none", se.fit = FALSE)
}
#' @export
.pred_Surv_lb <- function(y,w){
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
      class = 'survfit_lb'
    ))
  }

  U <- unique(y[,2][which(delta==1)])
  A <- y[,1]
  Z <- y[,2]
  M <- max(Z)
  V <- Z-A
  AV <- sort(unique(c(A,V)))

  dN <- sapply(U, function(s) sum(w*(Z==s)*delta))
  dQ <- sapply(AV, function(s) sum(w*(A==s)) + sum(w*(V==s)*delta))
  K <- sapply(AV, function(s) sum(w*(A>=s)) + sum(w*(V>=s)))
  K[K<=0] <- eps
  S_A <- sapply(U, function(s){
    prod(1-dQ[AV<=s]/K[AV<=s])
  })
  S_A[S_A<0] <- 0; S_A[S_A>1] <- 1
  R2 <- sapply(U, function(s) sum(w*(Z>=s))) - S_A*n
  R2[R2<=0] <- eps

  dN_R2 <- dN/R2
  dN_R2[dN_R2>1] <- 1

  S_pred <- sapply(Z, function(x){
    prod(1-dN_R2[U<=x])
  })

  S_pred[S_pred<0] <- 0; S_pred[S_pred>1] <- 1

  return(structure(
    list(time=Z, surv=S_pred),
    class = 'survfit_lb'
  ))
}
#' @export
print.survfit_lb <- function(x, ...){
  cat("LB Survival Object\n")
  cat("Median survival time :",format(median(x$time),digits=3))
}
