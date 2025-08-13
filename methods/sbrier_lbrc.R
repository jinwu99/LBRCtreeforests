# Implementation adapted from LTRCforests:
# https://github.com/weichiyao/TimeVaryingData_LTRCforests/blob/main/pkg/LTRCforests/R/sbrier_ltrc.R
# .ibsfunc and .bsfunc are unchanged
sbrier_lbrc <- function(obj, id = NULL, pred, type = c("IBS","BS")){
  if(!inherits(obj, "Surv"))
    stop("obj is not of class Surv")

  if (attr(obj, "type") != "counting")
    stop("only dataset with left-truncated (or length-biased) and right-censored (pseudo-subject) observations allowed")

  n <- nrow(obj)

  if (!is.null(id)){
    if (n != length(id)) stop("The length of id is different from the Surv object!")
  } else {
    id <- 1:n
  }

  id.sub = unique(id)
  n.sub = length(id.sub)

  obj <- as.data.frame(as.matrix(obj))

  if (n == n.sub){
    data_sbrier = obj
    data_sbrier$id = 1:n
  } else {
    data_sbrier <- data.frame(matrix(0, nrow = n.sub, ncol = 3))
    names(data_sbrier) <- c("start", "stop", "status")
    data_sbrier$id <- id.sub
    for (ii in 1:n.sub){
      data_sbrier[ii, ]$start = min(obj[id == id.sub[ii], ]$start)
      data_sbrier[ii, ]$stop = max(obj[id == id.sub[ii], ]$stop)
      data_sbrier[ii, ]$status = sum(obj[id == id.sub[ii], ]$status)
    }
  }

  #######################
  # Required for unbiased survival estimation under the no–left-truncation assumption
  data_sbrier$start <- 0
  #######################

  if (type[1] == "IBS"){
    ret <- sapply(1:n.sub, function(Ni) .ibsfunc(Ni = Ni, data_sbrier = data_sbrier, pred = pred))
    ret <- mean(ret)
    names(ret) = "Integrated Brier score"
  } else if (type[1] == "BS"){
    tpnt <- pred$survival.times[pred$survival.times <= min(pred$survival.tau)]
    bsres <- sapply(1:n.sub, function(Ni) .bsfunc(Ni = Ni, data_sbrier = data_sbrier, pred = pred, tpnt = tpnt))
    bsres <- rowMeans(bsres)
    ret <- data.frame(matrix(0, ncol = 2, nrow = length(tpnt)))
    colnames(ret) <- c("Time", "BScore")
    ret$Time <- tpnt
    ret$BScore <- bsres
  } else {
    stop("type can only be 'IBS' or 'BS'")
  }
  return(ret)
}

.ibsfunc <- function(Ni, data_sbrier, pred){
  id_uniq <- unique(data_sbrier$id)
  tpnt = pred$survival.times[pred$survival.times <= pred$survival.tau[Ni]]
  tlen = length(tpnt)
  ## Get the estimated survival probabilities
  if (class(pred$survival.probs)[1] == "matrix"){
    Shat = pred$survival.probs[1:tlen, Ni]
  } else if(class(pred$survival.probs)[1] == "list"){
    Shat = pred$survival.probs[[Ni]][1:tlen]
  }

  ######================ reverse Kaplan-Meier: estimate censoring distribution ====== ########
  # deal with ties
  hatcdist <- prodlim::prodlim(Surv(start, stop, status) ~ 1, data = data_sbrier, reverse = TRUE)

  Ttildei <- data_sbrier[data_sbrier$id == id_uniq[Ni], ]$stop

  Tleft = data_sbrier[data_sbrier$id == id_uniq[Ni], ]$start

  csurv_adj = predict(hatcdist, times = Tleft, type = "surv")
  if (is.na(csurv_adj)) stop("reverse Kaplan-Meier estimate at the left-truncateion point is NA! ")
  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type = "surv") / csurv_adj
  csurv_obs[csurv_adj == 0] <- Inf
  csurv_obs[csurv_obs == 0] <- Inf

  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv") / csurv_adj
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf

  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1/csurv_t, rep(1 / csurv_obs, sum(tpnt >= Ttildei)))

  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt < Ttildei)
  Indicator_t[Indicator_t == 0] = as.integer(data_sbrier[data_sbrier$id == id_uniq[Ni],]$status == 1)

  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt < Ttildei) - Shat) ^ 2 * csurv * Indicator_t
  ibs = diff(tpnt) %*% (fibs_itg[-length(fibs_itg)] + fibs_itg[-1]) / 2
  ibs = ibs / diff(range(tpnt))
  ibs
}

.bsfunc <- function(Ni, data_sbrier, pred, tpnt){
  id_uniq <- unique(data_sbrier$id)
  tlen = length(tpnt)
  ## Get the estimated survival probabilities
  if (class(pred$survival.probs)[1] == "matrix"){
    Shat = pred$survival.probs[1:tlen, Ni]
  } else if(class(pred$survival.probs)[1] == "list"){
    Shat = pred$survival.probs[[Ni]][1:tlen]
  }

  ######================ reverse Kaplan-Meier: estimate censoring distribution ================###########
  # deal with ties
  hatcdist <- prodlim::prodlim(Surv(start, stop, status) ~ 1, data = data_sbrier, reverse = TRUE)

  Ttildei <- data_sbrier[data_sbrier$id == id_uniq[Ni], ]$stop

  Tleft = data_sbrier[data_sbrier$id == id_uniq[Ni], ]$start

  csurv_adj = predict(hatcdist, times = Tleft, type  = "surv")
  if (is.na(csurv_adj)) stop("reverse Kaplan-Meier estimate at the left-truncateion point is NA! ")

  ### conditional survival for Observed value < t, G(Obs)
  csurv_obs <- predict(hatcdist, times = Ttildei, type  = "surv") / csurv_adj
  csurv_obs[csurv_adj == 0] <- Inf
  csurv_obs[csurv_obs == 0] <- Inf

  # conditional survival for Observed value > t, G(t)
  csurv_t <- predict(hatcdist, times = tpnt[tpnt < Ttildei], type = "surv") / csurv_adj
  csurv_t[is.na(csurv_t)] <- min(csurv_t, na.rm = TRUE)
  csurv_t[csurv_t == 0] <- Inf

  ## c(G^{-1}(t), G^{-1}(Obs))
  csurv <- c(1 / csurv_t, rep(1 / csurv_obs, sum(tpnt >= Ttildei)))

  ######================ indicator ================#################
  Indicator_t <- as.integer(tpnt < Ttildei)
  Indicator_t[Indicator_t == 0] = as.integer(data_sbrier[data_sbrier$id == id_uniq[Ni], ]$status == 1)

  ######================ Brier score =================#################
  fibs_itg = (as.integer(tpnt < Ttildei) - Shat) ^ 2 * csurv * Indicator_t
  fibs_itg
}
