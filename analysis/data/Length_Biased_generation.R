LBRC_sampling <- function(n, true_fail, true_params, ksi = 100, exp_cens_rate){
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

LBRC.generate_tree <- function(n=100, Dist = "WI", exp_cens_rate = 0.08, ksi = 100,
                               cov_set_num = 6, Test_mode = FALSE){
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # We will generate failure distribution conditional on specific covariate values
  Data$X1 <- sample(1:6, n, replace=T)
  Data$X2 <- sample(c(0,1), n, replace=T)
  Data$X3 <- runif(n,0,1)

  # False covariates
  for(cv in 1:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(1:6, n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n,0,1)
  }

  # true terminal node labels
  true_terminal_node_label <- c(3,4,6,7) # always
  Data[(Data$X1<=3 & Data$X2==0), 'tr_trmnd'] <- true_terminal_node_label[1]
  Data[(Data$X1<=3 & Data$X2==1),  'tr_trmnd'] <- true_terminal_node_label[2]
  Data[(Data$X1>3 & Data$X3<=0.5),  'tr_trmnd'] <- true_terminal_node_label[3]
  Data[(Data$X1>3 & Data$X3>0.5),   'tr_trmnd'] <- true_terminal_node_label[4]

  # term_nds : terminal nodes
  Dist_FUN_term_nds <- vector("list", 4)
  Dist_params_term_nds <- c()

  true_fail_dists <- function(Dist = "WI"){
    if(Dist == "Exp"){
      rdist <- rexp
      pdist <- pexp
      params_list <- list(
        list(rate = 0.1),
        list(rate = 0.23),
        list(rate = 0.4),
        list(rate = 0.9)
      )
    }else if(Dist == "WI"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(
        list(shape = 3, scale = 2.0),
        list(shape = 3, scale = 4.3),
        list(shape = 3, scale = 6.2),
        list(shape = 3, scale = 10.0)
      )
    }else if(Dist == "WD"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(
        list(shape = 0.9, scale = 7.0),
        list(shape = 0.9, scale = 3.0),
        list(shape = 0.9, scale = 2.5),
        list(shape = 0.9, scale = 1.0)
      )
    }else if(Dist == "Lgn"){
      pdist <- plnorm
      rdist <- rlnorm
      params_list <- list(
        list(meanlog = 2.0, sdlog = 0.3),
        list(meanlog = 1.7, sdlog = 0.2),
        list(meanlog = 1.3, sdlog = 0.3),
        list(meanlog = 0.5, sdlog = 0.5)
      )
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list))
  }

  Dist_info <- true_fail_dists(Dist = Dist)
  Dist_FUN_term_nds <- Dist_info$rdist
  Dist_params_term_nds <- Dist_info$params_list

  if(Test_mode == TRUE) exp_cens_rate <- 1e-100 # gives no censoring for evaluation

  # generate LBRC data for each terminal node
  true_fail <- Dist_FUN_term_nds
  for(nd in 1:4){
    true_params <- Dist_params_term_nds[[nd]]
    dat <- LBRC_sampling(n = sum(Data$tr_trmnd==true_terminal_node_label[nd]),
                         true_fail = true_fail,
                         true_params = true_params,
                         ksi = ksi,
                         exp_cens_rate = exp_cens_rate)

    Data[Data$tr_trmnd == true_terminal_node_label[nd],'A'] <- dat$left_trunc_time
    Data[Data$tr_trmnd == true_terminal_node_label[nd],'Y'] <- dat$event_time
    Data[Data$tr_trmnd == true_terminal_node_label[nd],'event'] <- dat$event
  }

  if(Test_mode == TRUE){
    ## time points of interest to evaluate the true survival function
    time.uniq <- unique(sort(c(0,Data$A,Data$Y,Data$Y-Data$A)))
    time.uniq <- time.uniq[time.uniq <= max(Data$Y)]
    tlen <- length(time.uniq)

    Dist_cdf <- Dist_info$pdist

    # true survival function given truncation time with no censoring
    true_surv_T_A <- matrix(NA,nrow=tlen, ncol=n)
    for(i in 1:n){
      true_surv_T_A_i <- vector(mode = "numeric", length=tlen)
      true_surv_T_A_i[time.uniq < Data$A[i]] <- 1
      # dist parameter of ith data
      ith_param <- Dist_params_term_nds[[match(Data$tr_trmnd[i],true_terminal_node_label)]]
      # P(T>t), t>=a_i
      S_T <- 1-do.call(Dist_cdf, c(list(q=time.uniq[time.uniq>=Data$A[i]]),ith_param))
      # P(T>a_i)
      S_A <- 1-do.call(Dist_cdf, c(list(q=Data$A[i]),ith_param))
      S_A <- max(S_A, 1e-10)
      true_surv_T_A_i[time.uniq >= Data$A[i]] <- S_T/S_A
      true_surv_T_A[,i] <- true_surv_T_A_i
    }
  }else{
    time.uniq <- NULL
    true_surv_T_A <- NULL
  }

  # delete true terminal node label from the Data
  tr_trmnd <- Data$tr_trmnd
  Data <- subset(Data, select = -tr_trmnd)

  # assign id
  Data$id <- 1:n

  cat(Dist, "tree failure distribution Generated ...\n")

  return(list(Data          = Data,
              tr_trmnd      = tr_trmnd,
              time.uniq     = time.uniq,
              true_surv_T_A = true_surv_T_A))
}


LBRC.generate_PH_nPH <- function(n=100, Dist = "WI", exp_cens_rate = 0.08, ksi = 100,
                                 cov_set_num = 6, nonlinear = FALSE, Test_mode = FALSE){
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # We will generate failure distribution conditional on specific covariate values
  Data$X1 <- sample(1:6, n, replace=T)
  Data$X2 <- sample(c(0,1), n, replace=T)
  Data$X3 <- runif(n,0,1)

  # False covariates
  for(cv in 1:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(1:6, n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n,0,1)
  }

  if(nonlinear == TRUE){
    loc_params <- -(cos(pi*(Data$X1+Data$X2)) + sqrt(Data$X1+Data$X2))
  }else{ # if not specified, linear
    loc_params <- -(Data$X1+Data$X2)
  }
  Params <- exp(loc_params)

  Dist_FUN_term_nds <- vector("list", 4)
  Dist_params_term_nds <- c()

  true_fail_dists <- function(Dist = "WI", Params){
    param_len <- length(Params)
    if(Dist == "Exp"){
      rdist <- rexp
      pdist <- pexp
      params_list <- lapply(1:param_len, function(l){
        list(rate = Params[l])
      })
    }else if(Dist == "WI"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 2, scale = 10*Params[l])
      })
    }else if(Dist == "WD"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 0.5, scale = 5*Params[l])
      })
    }else if(Dist == "Lgn"){
      pdist <- plnorm
      rdist <- rlnorm
      params_list <- lapply(1:param_len, function(l){
        list(meanlog = 1.5, sdlog = Params[l])
      })
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list))
  }

  Dist_info <- true_fail_dists(Dist = Dist, Params = Params)
  Dist_FUN_PH_nPH <- Dist_info$rdist
  Dist_params_PH_nPH <- Dist_info$params_list

  if(Test_mode == TRUE) exp_cens_rate <- 1e-10 # gives no censoring for evaluation

  # generate LBRC data
  true_fail <- Dist_FUN_PH_nPH
  dat <- t(sapply(1:n, function(i){
    unlist(
      LBRC_sampling(n = 1,
                    true_fail = true_fail,
                    true_params = Dist_params_PH_nPH[[i]],
                    ksi = ksi,
                    exp_cens_rate = exp_cens_rate)
    )
  }))
  Data$A <- dat[,"left_trunc_time"]
  Data$Y <- dat[,"event_time"]
  Data$event <- dat[,"event"]

  if(Test_mode == TRUE){
    ## time points of interest to evaluate the true survival function
    time.uniq <- unique(sort(c(0,Data$A,Data$Y,Data$Y-Data$A)))
    time.uniq <- time.uniq[time.uniq <= max(Data$Y)]
    tlen <- length(time.uniq)

    Dist_cdf <- Dist_info$pdist

    # true survival function given truncation time with no censoring
    true_surv_T_A <- matrix(NA,nrow=tlen, ncol=n)
    for(i in 1:n){
      true_surv_T_A_i <- vector(mode = "numeric", length=tlen)
      true_surv_T_A_i[time.uniq < Data$A[i]] <- 1
      # dist parameter of ith data
      ith_param <- Dist_params_PH_nPH[[i]]
      # P(T>t), t>=a_i
      S_T <- 1-do.call(Dist_cdf, c(list(q=time.uniq[time.uniq>=Data$A[i]]),ith_param))
      # P(T>a_i)
      S_A <- 1-do.call(Dist_cdf, c(list(q=Data$A[i]),ith_param))
      S_A <- max(S_A, 1e-10)
      true_surv_T_A_i[time.uniq >= Data$A[i]] <- S_T/S_A
      true_surv_T_A[,i] <- true_surv_T_A_i
    }
  }else{
    time.uniq <- NULL
    true_surv_T_A <- NULL
  }

  # assign id
  Data$id <- 1:n

  lab <- ifelse(nonlinear, "nonlinear", "linear")
  cat(Dist, lab, "failure distribution Generated ...\n")

  return(list(Data          = Data,
              time.uniq     = time.uniq,
              true_surv_T_A = true_surv_T_A))
}
