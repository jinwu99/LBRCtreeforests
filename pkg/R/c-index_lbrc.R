cindex_lbrc <- function(pred) {
  S <- pred$survival.probs
  risk <- colSums(1 - S, na.rm = TRUE)
  ok <- is.finite(risk)

  fit <- survival::concordance(
    pred$survival.obj[ok] ~ risk[ok],
    reverse = TRUE
  )

  fit$concordance
}
