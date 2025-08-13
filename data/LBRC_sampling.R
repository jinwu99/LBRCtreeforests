# Generate length-biased and right-censored (LBRC) data
LBRC_sampling <- function(n, true_fail, true_params, ksi = 500, exp_cens_rate){
  Count = 1
  ret = NULL

  while (Count <= n) {
    # truncation time
    w0 <- runif(1,0,ksi)
    # true failure time
    t0<- do.call(true_fail, c(list(n=1), true_params))
    if (t0 >= ksi-w0)
    {
      a<-ksi-w0
      t<-t0
      v<-t-a
      c<-rexp(1,exp_cens_rate)
      delta<-ifelse(v<c,1,0)
      y<-min(v,c)

      ret = rbind(ret, c(a, a+y, delta))
      Count <- Count + 1
    }
  }

  ret = data.frame(ret)
  dimnames(ret)[[2]] = c("left_trunc_time", "event_time","event")
  return(ret)
}
