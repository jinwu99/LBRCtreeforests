predict_onesample_cox <- function(pred, data, time.eval, target = "observed"){
  n <- nrow(data)
  time.tau <- rep(max(time.eval), n)
  if(target == "unbiased"){
    # Force truncation time to 0 to yield UNBIASEDL survival predictions
    data[, "A"] <- 0
  }
  pred$survival.probs <- sapply(1:n, function(i){
    .shatfunc(i, data = data, pred = pred, tpnt = time.eval, tau = time.tau)
  })
  return(pred)
}

Loss.func <- function(Est.Curves, True.Curves, T.pnt, time.tau, target = "observed"){
  if(target == "observed"){
    n <- ncol(Est.Curves)
    loss_i <- sapply(seq_len(n), function(i){
      tau_i <- time.tau[i]
      if(is.na(tau_i) || tau_i <= 0){
        return(NA_real_)
      }
      valid <- !is.na(Est.Curves[, i]) &
        !is.na(True.Curves[, i]) &
        T.pnt <= tau_i

      t_i <- T.pnt[valid]
      r_i <- (Est.Curves[valid, i] - True.Curves[valid, i])^2
      sum(diff(t_i) * (r_i[-1] + r_i[-length(r_i)]) / 2) / tau_i
    })

    return(mean(loss_i, na.rm = TRUE))
  }else{ # target == "unbiased"
      resid2_mat <- (Est.Curves-True.Curves)^2
      T.max <- max(T.pnt)
      T.num <- nrow(resid2_mat)
      n <- ncol(resid2_mat)
      return(sum( diff(T.pnt) * (resid2_mat[-1,] + resid2_mat[-T.num,])/2 ) / (T.max * n))
  }
}
