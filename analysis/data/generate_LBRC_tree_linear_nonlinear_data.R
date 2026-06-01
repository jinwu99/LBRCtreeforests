# Generate tree-structured LBRC data
LBRC.generate_tree <- function(n             = 100,
                               Dist          = "WI",
                               cens_rate     = 20,
                               ksi           = 500,
                               Test_mode     = FALSE,
                               cov_set_num   = 10,
                               target        = "observed"){
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # Failure distributions will be simulated based on these covariate values
  for(cv in 0:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(c(1:6), size=n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), size=n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n=n, min=0, max=1)
  }

  # true terminal node labels
  true_terminal_node_label <- c(1,2,3,4)
  Data[(Data$X1<=3 & Data$X2  ==0), 'tr_trmnd'] <- true_terminal_node_label[1]
  Data[(Data$X1<=3 & Data$X2  ==1), 'tr_trmnd'] <- true_terminal_node_label[2]
  Data[(Data$X1 >3 & Data$X3<=0.5), 'tr_trmnd'] <- true_terminal_node_label[3]
  Data[(Data$X1 >3 & Data$X3 >0.5), 'tr_trmnd'] <- true_terminal_node_label[4]

  # term_nds : terminal nodes
  Dist_FUN_term_nds <- vector("list", 4)
  Dist_params_term_nds <- c()

  true_fail_dists <- function(Dist = "WI", cens_rate = 20){
    if(Dist == "WI"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(
        list(shape = 2, scale = 10.0),
        list(shape = 2, scale = 6.0),
        list(shape = 2, scale = 3.5),
        list(shape = 2, scale = 2.0)
      )
      exp_cens_rates <- c(0.084, 0.34) # rate parameter of the exponential distribution for right-censoring
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }else if(Dist == "WD"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(
        list(shape = 0.9, scale = 7.0),
        list(shape = 0.9, scale = 3.0),
        list(shape = 0.9, scale = 2.5),
        list(shape = 0.9, scale = 1.0)
      )
      exp_cens_rates <- c(0.08, 0.36)
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }else if(Dist == "Lgn"){
      rdist <- rlnorm
      pdist <- plnorm
      params_list <- list(
        list(meanlog = 2.0, sdlog = 0.3),
        list(meanlog = 1.8, sdlog = 0.2),
        list(meanlog = 1.2, sdlog = 0.3),
        list(meanlog = 0.5, sdlog = 0.5)
      )
      exp_cens_rates <- c(0.092, 0.35)
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }else if(Dist == "Bat"){
      rdist <- rbathtub
      pdist <- pbathtub
      params_list <- list(
        list(a = 0.01),
        list(a = 0.06),
        list(a = 0.2),
        list(a = 0.7)
      )
      exp_cens_rates <- c(0.09, 0.39)
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list,
                exp_cens_rate = exp_cens_rate))
  }

  Dist_info <- true_fail_dists(Dist = Dist, cens_rate = cens_rate)
  Dist_FUN_term_nds <- Dist_info$rdist
  Dist_params_term_nds <- Dist_info$params_list
  exp_cens_rate <- Dist_info$exp_cens_rate

  if(Test_mode == TRUE & target == "unbiased"){
    exp_cens_rate <- 1e-100 # gives no censoring for evaluation
    ksi <- 0 # gives no truncation
  }

  # generate LBRC data for each terminal node
  while(1){
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
    if((length(unique(Data[,'Y'])) == n) & (min(Data[,'Y'])>0)) break
  }

  if(Test_mode == TRUE){
    if(target == "observed"){
      time.uniq <- c(0,sort(Data$Y))
      tlen <- length(time.uniq)
      Dist_cdf <- Dist_info$pdist

      # true survival function given truncation time
      true_surv <- matrix(NA,nrow=tlen, ncol=n)
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
        true_surv[,i] <- true_surv_T_A_i
      }
    }else{ # target == "unbiased"
      time.uniq <- c(0,sort((Data$Y)))
      tlen <- length(time.uniq)
      Dist_cdf <- Dist_info$pdist
      # true survival function
      true_surv <- matrix(NA,nrow=tlen, ncol=n)
      for(i in 1:n){
        # dist parameter of ith data
        ith_param <- Dist_params_term_nds[[match(Data$tr_trmnd[i],true_terminal_node_label)]]
        # P(T>t)
        S_T <- 1-do.call(Dist_cdf, c(list(q=time.uniq[time.uniq>=0]),ith_param))
        true_surv[,i] <- S_T
      }
    }
  }else{
    time.uniq <- NULL
    true_surv <- NULL
  }

  # delete true terminal node label from the Data
  tr_trmnd <- Data$tr_trmnd
  Data <- subset(Data, select = -tr_trmnd)

  # assign id
  Data$id <- 1:n

  return(list(Data          = Data,
              time.uniq     = time.uniq,
              true_surv     = true_surv))
}


# Generate (non)linear-structured LBRC data
LBRC.generate_PH_nPH <- function(n             = 100,
                                 Dist          = "WI",
                                 cens_rate     = 20,
                                 ksi           = 500,
                                 Test_mode     = FALSE,
                                 cov_set_num   = 10,
                                 model         = "linear",
                                 target        = "observed"){
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # Failure distributions will be simulated based on these covariate values
  for(cv in 0:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(c(1:6), size=n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), size=n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n=n, min=0, max=1)
  }

  if(model == "nonlinear"){
    beta <- c(1, 1, 1/6)
    if(Dist=="WI"){ # Weibull Increasing
      beta0 <- -log(10)
    }else{ # Weibull Decreasing (WD)
      beta0 <- -log(5)
    }
    loc_params <- (beta0 + beta[1]*cos(pi*(Data$X3+Data$X2)) +
                     beta[2]*sqrt(Data$X3+Data$X2) +
                     beta[3]*(Data$X1)^(Data$X2))
  }else if(model == "interaction"){
    beta <- c(0.5, -1.5)
    if(Dist=="WI"){
      beta0 <- -log(2)
      # loc_params <- (-log(2) + 1/2*Data$X1*Data$X2*Data$X3 - 1.5*Data$X3^3)
    }else{ # WD
      beta0 <- -log(1)
      # loc_params <- (-log(1) + 1/2*Data$X1*Data$X2*Data$X3 - 1.5*Data$X3^3)
    }
    loc_params <- (beta0 + beta[1]*Data$X1*Data$X2*Data$X3 + beta[2]*Data$X3^3)
  }else{ # if not specified, linear
    beta <- c(1, 1, -1/3)
    if(Dist=="WI"){ # Weibull Increasing
      beta0 <- -log(2);
    }else{ # Weibull Decreasing (WD)
      beta0 <- -log(1)
    }
    loc_params <- (beta0 + beta[1]*Data$X3 + beta[2]*Data$X2 + beta[3]*Data$X1)
  }

  Params <- exp(-loc_params)

  Dist_FUN_term_nds <- vector("list", 4)
  Dist_params_term_nds <- c()

  true_fail_dists <- function(Dist = "WI", Params, model, cens_rate = 20){
    param_len <- length(Params)
    if(Dist == "WI"){ # Weibull Increasing
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 2, scale = Params[l])
      })
      if(model == "nonlinear"){
        exp_cens_rates <- c(0.13, 0.6) # rate parameter of the exponential distribution for right-censoring
        if(cens_rate == 50){
          exp_cens_rate <- exp_cens_rates[2]
        }else{ # censoring rate 20%
          exp_cens_rate <- exp_cens_rates[1]
        }
      }else if(model == "interaction"){
        exp_cens_rates <- c(0.19, 0.8)
        if(cens_rate == 50){
          exp_cens_rate <- exp_cens_rates[2]
        }else{ # censoring rate 20%
          exp_cens_rate <- exp_cens_rates[1]
        }
      }else{ # linear model
        exp_cens_rates <- c(0.14, 0.65)
        if(cens_rate == 50){
          exp_cens_rate <- exp_cens_rates[2]
        }else{ # censoring rate 20%
          exp_cens_rate <- exp_cens_rates[1]
        }
      }
    }else{ # Weibull Decreasing (WD)
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 0.8, scale = Params[l])
      })
      if(model == "nonlinear"){
        exp_cens_rates <- c(0.13, 0.63)
        if(cens_rate == 50){
          exp_cens_rate <- exp_cens_rates[2]
        }else{ # censoring rate 20%
          exp_cens_rate <- exp_cens_rates[1]
        }
      }else if(model == "interaction"){
        exp_cens_rates <- c(0.17, 0.77)
        if(cens_rate == 50){
          exp_cens_rate <- exp_cens_rates[2]
        }else{ # censoring rate 20%
          exp_cens_rate <- exp_cens_rates[1]
        }
      }else{ # linear model
        exp_cens_rates <- c(0.12, 0.63)
        if(cens_rate == 50){
          exp_cens_rate <- exp_cens_rates[2]
        }else{ # censoring rate 20%
          exp_cens_rate <- exp_cens_rates[1]
        }
      }
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list,
                exp_cens_rate = exp_cens_rate))
  }

  Dist_info <- true_fail_dists(Dist = Dist, Params = Params, model = model, cens_rate = cens_rate)
  Dist_FUN_PH_nPH <- Dist_info$rdist
  Dist_params_PH_nPH <- Dist_info$params_list
  exp_cens_rate <- Dist_info$exp_cens_rate

    if(Test_mode == TRUE & target == "unbiased"){
    exp_cens_rate <- 1e-100 # gives no censoring for evaluation
    ksi <- 0 # gives no truncation
  }

  while(1){ # keep generating LBRC data until unique time points match the sample size
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
    if((length(unique(Data[,'Y']))==n) & (min(Data[,'Y'])>0)) break
  }

  if(Test_mode == TRUE){
    if(target == "observed"){
      time.uniq <- c(0,sort(Data$Y))
      tlen <- length(time.uniq)
      Dist_cdf <- Dist_info$pdist

      # true survival function given truncation time
      true_surv <- matrix(NA,nrow=tlen, ncol=n)
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
        true_surv[,i] <- true_surv_T_A_i
      }
    }else{ # target == "unbiased"
      ## time points of interest to evaluate the true survival function
      time.uniq <- c(0,sort((Data$Y)))
      tlen <- length(time.uniq)
      Dist_cdf <- Dist_info$pdist
      # true survival function
      true_surv <- matrix(NA,nrow=tlen, ncol=n)
      for(i in 1:n){
        # dist parameter of ith data
        ith_param <- Dist_params_PH_nPH[[i]]
        # P(T>t)
        S_T <- 1-do.call(Dist_cdf, c(list(q=time.uniq[time.uniq>=0]),ith_param))
        true_surv[,i] <- S_T
      }
    }
  }else{
    time.uniq <- NULL
    true_surv <- NULL
  }

  # assign id
  Data$id <- 1:n

  cat(Dist, model, "failure dist with N",n,"P",cov_num,
      "C",mean(1-Data$event),"Generated ...\n")

  return(list(Data          = Data,
              time.uniq     = time.uniq,
              true_surv     = true_surv))
}


# Sensitivity analysis: generate general LTRC tree or nonlinear data
LTRC.generate_tree_nPH <- function(n             = 200,
                                   rho           = 0,
                                   tau           = 40,
                                   Test_mode     = FALSE,
                                   cov_set_num   = 10,
                                   scenario      = "tree_texpt"){ # tree_covd, nlin_texpt, nlin_covd
  trunc_rate <- round(rho/tau, 2)

  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # Failure distributions will be simulated based on these covariate values
  for(cv in 0:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(c(1:6), size=n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), size=n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n=n, min=0, max=1)
  }

  if(scenario %in% c("tree_texpt", "tree_covd")){
    # true terminal node labels
    true_terminal_node_label <- c(1,2,3,4)
    Data[(Data$X1<=3 & Data$X2  ==0), 'tr_trmnd'] <- true_terminal_node_label[1]
    Data[(Data$X1<=3 & Data$X2  ==1), 'tr_trmnd'] <- true_terminal_node_label[2]
    Data[(Data$X1 >3 & Data$X3<=0.5), 'tr_trmnd'] <- true_terminal_node_label[3]
    Data[(Data$X1 >3 & Data$X3 >0.5), 'tr_trmnd'] <- true_terminal_node_label[4]

    # term_nds : terminal nodes
    Dist_FUN_term_nds <- vector("list", 4)
    Dist_params_term_nds <- c()

    true_fail_dists <- function(scenario, trunc_rate){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(
        list(shape = 2, scale = 10.0),
        list(shape = 2, scale = 6.0),
        list(shape = 2, scale = 3.5),
        list(shape = 2, scale = 2.0)
      )
      if(scenario == "tree_texpt"){
        exp_cens_rate <- switch(as.character(trunc_rate),
                                "0.1"=0.079, "0.2"=0.077, "0.5"=0.07, "1"=0.06)
      }else{ # scenario == "tree_covd"
        exp_cens_rate <- switch(as.character(trunc_rate),
                                "0.1"=0.079, "0.2"=0.077, "0.5"=0.07, "1"=0.06)
      }
      return(list(rdist = rdist,
                  pdist = pdist,
                  params_list = params_list,
                  exp_cens_rate = exp_cens_rate))
    }

    Dist_info <- true_fail_dists(scenario = scenario, trunc_rate = trunc_rate)
    Dist_FUN_term_nds <- Dist_info$rdist
    Dist_params_term_nds <- Dist_info$params_list
    exp_cens_rate <- Dist_info$exp_cens_rate

    # generate LTRC data for each terminal node
    if(scenario == "tree_covd"){
      while(1){
        true_fail <- Dist_FUN_term_nds
        dat <- t(sapply(1:n, function(i){
          unlist(
            LTRC_sampling(n = 1,
                          true_fail = true_fail,
                          true_params = Dist_params_term_nds[[Data[i,"tr_trmnd"]]],
                          exp_cens_rate = exp_cens_rate,
                          rho = ifelse(Data[i,"X1"]<=3, rho, 0), tau = tau)
          )
        }))
        Data$A <- dat[,"left_trunc_time"]
        Data$Y <- dat[,"event_time"]
        Data$event <- dat[,"event"]
        if((length(unique(Data[,'Y']))==n) & (min(Data[,'Y'])>0)) break
      }
    }else{ # tree_texpt
      while(1){
        true_fail <- Dist_FUN_term_nds
        for(nd in 1:4){
          true_params <- Dist_params_term_nds[[nd]]
          dat <- LTRC_sampling(n = sum(Data$tr_trmnd==true_terminal_node_label[nd]),
                               true_fail = true_fail,
                               true_params = true_params,
                               exp_cens_rate = exp_cens_rate,
                               rho = rho, tau = tau)
          Data[Data$tr_trmnd == true_terminal_node_label[nd],'A'] <- dat$left_trunc_time
          Data[Data$tr_trmnd == true_terminal_node_label[nd],'Y'] <- dat$event_time
          Data[Data$tr_trmnd == true_terminal_node_label[nd],'event'] <- dat$event
        }
        if((length(unique(Data[,'Y'])) == n) & (min(Data[,'Y'])>0)) break
      }
    }

    if(Test_mode == TRUE){
      ## time points of interest to evaluate the true survival function
      # tlen <- 1000
      # time.uniq <- seq(0,max(Data$Y),length.out=tlen)
      time.uniq <- c(0,sort(Data$Y))
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
  }else if(scenario %in% c("nlin_texpt", "nlin_covd")){
    beta <- c(1, 1, 1/6)
    beta0 <- -log(10)
    loc_params <- (beta0 + beta[1]*cos(pi*(Data$X3+Data$X2)) +
                     beta[2]*sqrt(Data$X3+Data$X2) +
                     beta[3]*(Data$X1)^(Data$X2))
    Params <- exp(-loc_params)

    Dist_FUN_term_nds <- vector("list", 4)
    Dist_params_term_nds <- c()

    true_fail_dists <- function(scenario, trunc_rate){
      param_len <- length(Params)
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 2, scale = Params[l])
      })
      if(scenario == "nlin_texpt"){
        exp_cens_rate <- switch(as.character(trunc_rate),
                                "0.1"=0.11, "0.2"=0.11, "0.5"=0.10, "1"=0.09)
      }else{ # scenario == "nlin_covd"
        exp_cens_rate <- switch(as.character(trunc_rate),
                                "0.1"=0.12, "0.2"=0.12, "0.5"=0.11, "1"=0.10)
      }
      return(list(rdist = rdist,
                  pdist = pdist,
                  params_list = params_list,
                  exp_cens_rate = exp_cens_rate))
    }

    Dist_info <- true_fail_dists(scenario = scenario, trunc_rate = trunc_rate)
    Dist_FUN_PH_nPH <- Dist_info$rdist
    Dist_params_PH_nPH <- Dist_info$params_list
    exp_cens_rate <- Dist_info$exp_cens_rate

    if(scenario == "nlin_texpt"){
      while(1){ # keep generating LTRC data until unique time points match the sample size
        true_fail <- Dist_FUN_PH_nPH
        dat <- t(sapply(1:n, function(i){
          unlist(
            LTRC_sampling(n = 1,
                          true_fail = true_fail,
                          true_params = Dist_params_PH_nPH[[i]],
                          exp_cens_rate = exp_cens_rate,
                          rho = rho, tau = tau)
          )
        }))
        Data$A <- dat[,"left_trunc_time"]
        Data$Y <- dat[,"event_time"]
        Data$event <- dat[,"event"]
        if((length(unique(Data[,'Y']))==n) & (min(Data[,'Y'])>0)) break
      }
    }else{ # scenario == "nlin_covd"
      while(1){ # keep generating LTRC data until unique time points match the sample size
        true_fail <- Dist_FUN_PH_nPH
        dat <- t(sapply(1:n, function(i){
          unlist(
            LTRC_sampling(n = 1,
                          true_fail = true_fail,
                          true_params = Dist_params_PH_nPH[[i]],
                          exp_cens_rate = exp_cens_rate,
                          rho = ifelse(Data[i,"X1"]<=3, rho, 0), tau = tau)
          )
        }))
        Data$A <- dat[,"left_trunc_time"]
        Data$Y <- dat[,"event_time"]
        Data$event <- dat[,"event"]
        if((length(unique(Data[,'Y']))==n) & (min(Data[,'Y'])>0)) break
      }
    }

    if(Test_mode == TRUE){
      ## time points of interest to evaluate the true survival function
      # tlen <- 1000
      # time.uniq <- seq(0,max(Data$Y),length.out=tlen)
      time.uniq <- c(0,sort(Data$Y))
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
  }

  # assign id
  Data$id <- 1:n

  lab <- scenario
  cat("WI", lab, "failure dist with N",n,"P",cov_num,
      "C",mean(1-Data$event),"Generated ...\n")

  return(list(Data          = Data,
              time.uniq     = time.uniq,
              true_surv     = true_surv_T_A))
}



# Generate one-sample LBRC data independent of covariates
# This is for testing of unbiasedness of LBRC-conditional inference tree
LBRC.generate_one_sample <- function(n         = 200,
                                     Dist      = "WI",
                                     cens_rate = 20,
                                     ksi       = 500){
  cov_set_num <- 2
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate fake covariate values in advance
  for(cv in 0:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(c(1:6), size=n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), size=n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n=n, min=0, max=1)
  }


  true_fail_dists <- function(Dist = "WI", cens_rate = 20){
    if(Dist == "WI"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(shape = 2, scale = 3.0)
      exp_cens_rates <- c(0.15, 0.5) # rate parameter of the exponential distribution for right-censoring
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }else if(Dist == "WD"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(shape = 0.9, scale = 2.0)
      exp_cens_rates <- c(0.11, 0.44)
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }else if(Dist == "Lgn"){
      rdist <- rlnorm
      pdist <- plnorm
      params_list <- list(meanlog = 1.4, sdlog = 0.4)
      exp_cens_rates <- c(0.09, 0.33)
      if(cens_rate == 50){
        exp_cens_rate <- exp_cens_rates[2]
      }else{ # censoring rate 20%
        exp_cens_rate <- exp_cens_rates[1]
      }
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list,
                exp_cens_rate = exp_cens_rate))
  }

  Dist_info <- true_fail_dists(Dist = Dist, cens_rate = cens_rate)

  dat <- LBRC_sampling(n = n,
                       true_fail = Dist_info$rdist,
                       true_params = Dist_info$params_list,
                       ksi = ksi,
                       exp_cens_rate = Dist_info$exp_cens_rate)

  Data$A <- dat[,"left_trunc_time"]
  Data$Y <- dat[,"event_time"]
  Data$event <- dat[,"event"]

  Data$id <- 1:n

  return(Data)
}



LTRC.generate_one_sample <- function(n         = 200,
                                     rho       = 0,
                                     tau       = 40,
                                     scenario  = "unbias_texpt"){
  trunc_rate <- round(rho/tau, 2)

  cov_set_num <- 2
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate fake covariate values in advance
  for(cv in 0:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- sample(c(1:6), size=n, replace=T)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), size=n, replace=T)
    Data[,paste0("X",3*cv+3)] <- runif(n=n, min=0, max=1)
  }

  true_fail_dists <- function(scenario, trunc_rate){
    rdist <- rweibull
    pdist <- pweibull
    params_list <- list(shape = 2, scale = 3.0)
    if(scenario == "unbias_texpt"){
      exp_cens_rate <- switch(as.character(trunc_rate),
                              "0.1"=0.13, "0.2"=0.13, "0.5"=0.12, "1"=0.11)
    }else{ # scenario == "unbias_covd"
      exp_cens_rate <- switch(as.character(trunc_rate),
                              "0.1"=0.14, "0.2"=0.14, "0.5"=0.13, "1"=0.12)
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list,
                exp_cens_rate = exp_cens_rate))
  }

  Dist_info <- true_fail_dists(scenario = scenario, trunc_rate = trunc_rate)

  if(scenario == "unbias_texpt"){
    dat <- LTRC_sampling(n = n,
                         true_fail = Dist_info$rdist,
                         true_params = Dist_info$params_list,
                         exp_cens_rate = Dist_info$exp_cens_rate,
                         rho = rho, tau = tau)
  }else{ # scenario == "unbias_covd"
    dat <- t(sapply(1:n, function(i){
      unlist(
        LTRC_sampling(n = 1,
                      true_fail = Dist_info$rdist,
                      true_params = Dist_info$params_list,
                      exp_cens_rate = Dist_info$exp_cens_rate,
                      rho = ifelse(Data[i,"X1"]<=3, rho, 0), tau = tau)
      )
    }))
  }

  Data$A <- dat[,"left_trunc_time"]
  Data$Y <- dat[,"event_time"]
  Data$event <- dat[,"event"]

  Data$id <- 1:n

  return(Data)
}






