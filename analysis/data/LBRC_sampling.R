# Generate length-biased and right-censored (LBRC) data
# Sampling scheme adopted from Paper:
#   Nonparametric estimation for length-biased and right-censored data (Huang and Qin, 2011)
LBRC_sampling <- function(n, true_fail, true_params, ksi = 500, exp_cens_rate){
  Count <- 1
  ret <- NULL

  while (Count <= n) {
    # onset time
    w0 <- runif(1, 0, ksi)
    # true failure time from onset time
    t0 <- do.call(true_fail, c(list(n = 1), true_params))
    if (t0 >= ksi - w0)
    {
      a <- ksi - w0
      t <- t0
      v <- t - a
      c <- rexp(1, exp_cens_rate)
      delta <- ifelse(v < c, 1, 0)
      y <- min(v, c)

      ret <- rbind(ret, c(a, a+y, delta))
      Count <- Count + 1
    }
  }

  ret = data.frame(ret)
  dimnames(ret)[[2]] = c("left_trunc_time", "event_time", "event")
  return(ret)
}

# Generate general left-truncated and right-censored (LTRC) data - for sensitivity analysis
LTRC_sampling <- function(n, true_fail, true_params, exp_cens_rate,
                          rho, tau){
  Count <- 1
  ret <- NULL

  while (Count <= n) {
    # true truncation time (assuming onset happened at time zero)
    a0 <- rtexp_tilt(1, rho, tau)
    # true failure time
    t0 <- do.call(true_fail, c(list(n = 1), true_params))
    if (t0 >= a0)
    {
      a <- a0
      t <- t0
      v <- t - a
      c <- rexp(1, exp_cens_rate)
      delta <- ifelse(v < c, 1, 0)
      y <- min(v, c)

      ret <- rbind(ret, c(a, a+y, delta))
      Count <- Count + 1
    }
  }

  ret = data.frame(ret)
  dimnames(ret)[[2]] = c("left_trunc_time", "event_time", "event")
  return(ret)
}
