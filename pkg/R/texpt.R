rtexp_tilt <- function(n, rho, tau) {
  if (rho == 0) {
    return(runif(n, 0, tau))  # uniform limit
  }
  u <- runif(n)
  a <- -tau / rho * log(1 - u * (1 - exp(-rho)))
  return(a)
}

dtexp_tilt <- function(x, rho, tau, log = FALSE) {
  dens <- numeric(length(x))

  inside <- (x >= 0) & (x <= tau)

  if (rho == 0) {
    dens[inside] <- 1 / tau
  } else {
    dens[inside] <- (rho / tau) * exp(-rho * x[inside] / tau) / (1 - exp(-rho))
  }

  if (log) {
    dens[inside] <- log(dens[inside])
    dens[!inside] <- -Inf
  } else {
    dens[!inside] <- 0
  }

  return(dens)
}


