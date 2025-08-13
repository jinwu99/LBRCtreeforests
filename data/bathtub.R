# from https://github.com/weichiyao/ICdata_ICcforest/blob/master/analysis/data/bathtub.R

###-----Inverse transform sampling-------------------######

endsign <- function(f, sign = 1) {
  b <- sign
  while (sign * f(b) < 0) b <- 10 * b
  return(b)
}

##=======================================================
rbathtub <- function(n, a, spdf.lower = 0, spdf.upper = Inf) {
  cdf <- function(x){1-exp(-a*x^2/2)/(1+5*x)^(1/5)}

  invcdf <- function(u){
    subcdf <- function(t) cdf(t) - u
    if (spdf.lower == -Inf)
      spdf.lower <- endsign(subcdf, -1)
    if (spdf.upper == Inf)
      spdf.upper <- endsign(subcdf)
    return(uniroot(subcdf, c(spdf.lower, spdf.upper))$root)
  }

  sapply(runif(n), invcdf)
}

##====================================================
pbathtub <- function(q, a) {
  1-exp(-a*q^2/2)/(1+5*q)^(1/5)
}


