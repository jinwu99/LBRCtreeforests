
predict_ctree_cox <- function(pred, data, time.eval){
  time.tau <- rep(max(time.eval), n)
  pred$survival.probs <- sapply(1:n, function(i){
    .shatfunc(i, data = data, pred = pred, tpnt = time.eval, tau = time.tau)
  })
  return(pred)
}

Loss.func <- function(Est.Curves, True.Curves, T.pnt){
  # dimension of Est.Curves and True.Curves :
  # length(T.pnt) x nunmber of data
  # row for each column(data) : 0 to maximum of T.pnt
  resid2_mat <- (Est.Curves-True.Curves)^2

  # integrate
  T.max <- max(T.pnt)
  T.num <- nrow(resid2_mat)
  n <- ncol(resid2_mat)

  return(sum( diff(T.pnt) * (resid2_mat[-1,] + resid2_mat[-T.num,])/2 ) / (T.max * n))
  # return(sum( (resid2_mat[-1,] + resid2_mat[-T.num,])/2 ) / (T.max * n))
}
