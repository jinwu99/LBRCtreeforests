# Plot method for constparty with LTRC/LBRC extensions
# Extends partykit::plot.constparty to support LTRC–CIT and LBRC–CIT variants
plot.constparty <- function(x, main = NULL,
                            terminal_panel = NULL, tp_args = list(),
                            inner_panel = node_inner, ip_args = list(),
                            edge_panel = edge_simple, ep_args = list(),
                            type = c("extended", "simple"), drop_terminal = NULL, tnex = NULL,
                            newpage = TRUE, pop = TRUE, gp = gpar(), ...)
{
  type <- match.arg(type)
  if (type == "simple") {
    x <- as.simpleparty(x)
    if (is.null(terminal_panel))
      terminal_panel <- node_terminal
    if (is.null(tnex)) tnex <- 1
    if (is.null(drop_terminal)) drop_terminal <- FALSE
    if (is.null(tp_args) || length(tp_args) < 1L) {
      tp_args <- list(FUN = .make_formatinfo_simpleparty(x, digits = getOption("digits") - 4L, sep = "\n"))
    } else {
      if(is.null(tp_args$FUN)) {
        tp_args$FUN <- .make_formatinfo_simpleparty(x, digits = getOption("digits") - 4L, sep = "\n")
      }
    }
  } else {
    if (is.null(terminal_panel)) {
      cl <- class(x$fitted[["(response)"]])
      if("factor" %in% cl) {
        terminal_panel <- node_barplot
      } else if("Surv" %in% cl) {
        surv_type <- x$perm_test_est
        if(surv_type == "MFLE"){
          main <- "LBRC-CIT-F"
        }else if(surv_type == "MCLE"){
          main <- "LBRC-CIT-C"
        }else{ # surv_type == "KM"
          main <- "LTRC-CIT"
        }
        terminal_panel <- node_surv2
      } else if ("data.frame" %in% cl) {
        terminal_panel <- node_mvar
        if (is.null(tnex)) tnex <- 2 * NCOL(x$fitted[["(response)"]])
      } else {
        terminal_panel <- node_boxplot
      }
    }
    if (is.null(tnex)) tnex <- 2
    if (is.null(drop_terminal)) drop_terminal <- TRUE
  }

  plot.party(x, main = main,
             terminal_panel = terminal_panel, tp_args = tp_args,
             inner_panel = inner_panel, ip_args = ip_args,
             edge_panel = edge_panel, ep_args = ep_args,
             drop_terminal = drop_terminal, tnex = tnex,
             newpage = newpage, pop = pop, gp = gp, ...)
}


node_surv2 <- function (obj, col = "black", bg = "white", yscale = c(0, 1),
                        ylines = 2, id = TRUE, mainlab = NULL, gp = gpar(), ...)
{
  y <- obj$fitted[["(response)"]]
  stopifnot(inherits(y, "Surv"))
  kmsurvfit <- function(y, weights, ...) survfit(y ~ 1, weights = weights)
  MCLE <- function(y, w){
    if (length(y) == 0) return(NA)
    idx = which(w>0)
    y = y[idx,]
    w = w[idx]

    n <- sum(w)
    delta <- y[,3]
    Z <- y[,2]

    if(sum(delta)==0){
      S_pred <- as.double(rep(1,dim(y)[1]))
      RET <- list(time=as.numeric(Z),
                  surv=as.numeric(S_pred))
      return(RET)
    }

    Z_order <- order(y[,2])
    delta <- y[Z_order,3]
    U <- unique(y[Z_order,2][which(delta==1)])
    A <- y[Z_order,1]
    Z <- y[Z_order,2]
    V <- Z-A

    dN <- sapply(U, function(s) sum(w*(Z==s)*delta))
    R <- sapply(U, function(s) sum( w*( (A<=s & Z>=s) + delta*(V<=s & Z>=s) ) ) )/2

    dN_R <- dN/R
    dN_R[dN_R>1] <- 1

    S_pred <- sapply(Z, function(x){
      prod(1-dN_R2[U<=x])
    })

    S_pred[S_pred<0] <- 0; S_pred[S_pred>1] <- 1

    RET <- list(time=as.numeric(Z),
                surv=as.numeric(S_pred))

    return(RET)
  }
  MFLE <- function(y,w,eps=1e-7,max_iter=100){
    if (length(y) == 0) return(NA)
    idx = which(w>0)
    y = y[idx,]
    w = w[idx]

    n <- sum(w)
    delta <- y[,3]
    Z <- y[,2]

    if(sum(delta)==0){
      S_pred <- as.double(rep(1,dim(y)[1]))
      RET <- list(time=as.numeric(Z),
                  surv=as.numeric(S_pred))
      return(RET)
    }

    res <- vardiCpp(y,w,eps = eps,max_iter = max_iter)
    Y <- res$t
    S_pred <- res$S
    S_pred[S_pred<0] <- 0; S_pred[S_pred>1] <- 1
    RET <- list(time=as.numeric(Y),
                surv=as.numeric(S_pred))

    return(RET)
  }
  dostep <- function(x, y) {
    if (is.na(x[1] + y[1])) {
      x <- x[-1]
      y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {
      dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
      n2 <- sum(dupy)
      xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
      yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
      RET <- list(x = xrep, y = yrep)
    }
    else {
      if (n == 1) {
        RET <- list(x = x, y = y)
      }
      else {
        RET <- list(x = x[c(1, 2, 2)], y = y[c(1, 1,
                                               2)])
      }
    }
    return(RET)
  }
  rval <- function(node) {
    nid <- id_node(node)
    dat <- data_party(obj, nid)
    yn <- dat[["(response)"]]
    wn <- dat[["(weights)"]]
    if (is.null(wn))
      wn <- rep(1, NROW(yn))
    surv_type <- obj$perm_test_est
    if (surv_type == "MFLE"){
      fit <- MFLE(yn, wn)
    }else if (surv_type == "MCLE"){
      fit <- MCLE(yn, wn)
    }else{
      fit <- kmsurvfit(yn, weights = wn, ...)
    }
    a <- dostep(fit$time, fit$surv)
    yscale <- yscale
    xscale <- c(0, max(y[, 1], na.rm = TRUE))
    top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                                            widths = unit(c(ylines, 1, 1), c("lines", "null",
                                                                             "lines")), heights = unit(c(1, 1), c("lines",
                                                                                                                  "null"))), width = unit(1, "npc"), height = unit(1,
                                                                                                                                                                   "npc") - unit(2, "lines"), name = paste("node_surv",
                                                                                                                                                                                                           nid, sep = ""), gp = gp)
    pushViewport(top_vp)
    grid.rect(gp = gpar(fill = bg, col = 0))
    top <- viewport(layout.pos.col = 2, layout.pos.row = 1)
    pushViewport(top)
    if (is.null(mainlab)) {
      mainlab <- if (id) {
        function(id, nobs) sprintf("Node %s (n = %s)",
                                   id, nobs)
      }
      else {
        function(id, nobs) sprintf("n = %s", nobs)
      }
    }
    if (is.function(mainlab)) {
      mainlab <- mainlab(nid, sum(wn))
    }
    grid.text(mainlab)
    popViewport()
    plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                     xscale = xscale, yscale = yscale, name = paste0("node_surv",
                                                                     nid, "plot"), clip = FALSE)
    pushViewport(plot)
    grid.xaxis()
    grid.yaxis()
    grid.rect(gp = gpar(fill = "transparent"))
    grid.clip()
    grid.lines(unit(a$x, "native"), unit(a$y, "native"),
               gp = gpar(col = col))
    upViewport(2)
  }
  return(rval)
}
class(node_surv2) <- "grapcon_generator"
