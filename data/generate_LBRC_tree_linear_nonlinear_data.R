# Generate tree-structured LBRC data
LBRC.generate_tree <- function(n             = 100,
                               Dist          = "WI",
                               exp_cens_rate = 1e-100,
                               ksi           = 500,
                               Test_mode     = FALSE,
                               cov_set_num   = 10){
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # We will generate failure distribution conditional on specific covariate values
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

  true_fail_dists <- function(Dist = "WI"){
    if(Dist == "WI"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(
        list(shape = 2, scale = 10.0),
        list(shape = 2, scale = 6.0),
        list(shape = 2, scale = 3.5),
        list(shape = 2, scale = 2.0)
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
      rdist <- rlnorm
      pdist <- plnorm
      params_list <- list(
        list(meanlog = 2.0, sdlog = 0.3),
        list(meanlog = 1.8, sdlog = 0.2),
        list(meanlog = 1.2, sdlog = 0.3),
        list(meanlog = 0.5, sdlog = 0.5)
      )
    }else if(Dist == "Bat"){
      rdist <- rbathtub
      pdist <- pbathtub
      params_list <- list(
        list(a = 0.01),
        list(a = 0.06),
        list(a = 0.2),
        list(a = 0.7)
      )
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list))
  }

  Dist_info <- true_fail_dists(Dist = Dist)
  Dist_FUN_term_nds <- Dist_info$rdist
  Dist_params_term_nds <- Dist_info$params_list

  if(Test_mode == TRUE){
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
    ## time points of interest to evaluate the true survival function
    time.uniq <- c(0,sort((Data$Y)))
    tlen <- length(time.uniq)
    Dist_cdf <- Dist_info$pdist
    # true survival function given truncation time with no censoring
    true_surv <- matrix(NA,nrow=tlen, ncol=n)
    for(i in 1:n){
      # dist parameter of ith data
      ith_param <- Dist_params_term_nds[[match(Data$tr_trmnd[i],true_terminal_node_label)]]
      # P(T>t)
      S_T <- 1-do.call(Dist_cdf, c(list(q=time.uniq[time.uniq>=0]),ith_param))
      true_surv[,i] <- S_T
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
                                 exp_cens_rate = 1e-100,
                                 ksi           = 500,
                                 Test_mode     = FALSE,
                                 cov_set_num   = 7,
                                 nonlinear     = FALSE){
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # Generate covariate values in advance
  # We will generate failure distribution conditional on specific covariate values
  for(cv in 0:(cov_set_num-1)){
    Data[,paste0("X",3*cv+1)] <- runif(n=n, min=0, max=1)
    Data[,paste0("X",3*cv+2)] <- sample(c(0,1), size=n, replace=T)
    Data[,paste0("X",3*cv+3)] <- sample(c(1:6), size=n, replace=T)
  }

  if(nonlinear == TRUE){
    beta <- c(1, 1, 1/6)
    if(Dist=="WI"){ # Weibull Increasing
      beta0 <- -log(10)
    }else{ # Weibull Decreasing (WD)
      beta0 <- -log(5)
    }
    loc_params <- (beta0 + beta[1]*cos(pi*(Data$X1+Data$X2)) +
                     beta[2]*sqrt(Data$X1+Data$X2) +
                     beta[3]*(Data$X3)^(Data$X2))
  }else{ # if not specified, linear
    beta <- c(1, 1, -1/3)
    if(Dist=="WI"){ # Weibull Increasing
      beta0 <- -log(2);
    }else{ # Weibull Decreasing (WD)
      beta0 <- -log(1)
    }
    loc_params <- (beta0 + beta[1]*Data$X1 + beta[2]*Data$X2 + beta[3]*Data$X3)
  }
  Params <- exp(-loc_params)

  Dist_FUN_term_nds <- vector("list", 4)
  Dist_params_term_nds <- c()

  true_fail_dists <- function(Dist = "WI", Params){
    param_len <- length(Params)
    if(Dist == "WI"){ # Weibull Increasing
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 2, scale = Params[l])
      })
    }else{ # Weibull Decreasing (WD)
      rdist <- rweibull
      pdist <- pweibull
      params_list <- lapply(1:param_len, function(l){
        list(shape = 0.8, scale = Params[l])
      })
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list))
  }

  Dist_info <- true_fail_dists(Dist = Dist, Params = Params)
  Dist_FUN_PH_nPH <- Dist_info$rdist
  Dist_params_PH_nPH <- Dist_info$params_list

  if(Test_mode == TRUE){
    exp_cens_rate <- 1e-100 # gives no censoring for evaluation
    ksi <- 0 # gives no truncation
  }

  while(1){ # generate LBRC data until the length of unique timepoints are same as sample size
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
    ## time points of interest to evaluate the true survival function
    time.uniq <- c(0,sort((Data$Y)))
    tlen <- length(time.uniq)
    Dist_cdf <- Dist_info$pdist
    # true survival function given truncation time with no censoring
    true_surv <- matrix(NA,nrow=tlen, ncol=n)
    for(i in 1:n){
      # dist parameter of ith data
      ith_param <- Dist_params_PH_nPH[[i]]
      # P(T>t)
      S_T <- 1-do.call(Dist_cdf, c(list(q=time.uniq[time.uniq>=0]),ith_param))
      true_surv[,i] <- S_T
    }
  }else{
    time.uniq <- NULL
    true_surv <- NULL
  }

  # assign id
  Data$id <- 1:n

  lab <- ifelse(nonlinear, "nonlinear", "linear")
  cat(Dist, lab, "failure dist with N",n,"P",cov_num,
      "C",mean(1-Data$event),"Generated ...\n")

  return(list(Data          = Data,
              time.uniq     = time.uniq,
              true_surv     = true_surv))
}


# Generate one-sample LBRC data - no association to any covariates
LBRC.generate_one_sample <- function(Dist="WD", n=300, exp_cens_rate=0.08, ksi=100){
  cov_set_num <- 2
  cov_num <- 3*cov_set_num
  col_num <- cov_num + 3 # A, Y, event
  Data <- as.data.frame(matrix(NA,n,col_num))
  names(Data) <- c(paste0("X",1:cov_num),'A','Y','event')

  # fake covariates
  Data[,"X1"] <- runif(n=n, min=0, max=1)
  Data[,"X2"] <- runif(n=n, min=0, max=1)
  Data[,"X3"] <- sample(seq(from=0,to=1,length.out=11),size=n,replace=T)
  Data[,"X4"] <- sample(seq(from=0,to=1,length.out=11),size=n,replace=T)
  Data[,"X5"] <- sample(c(0,1),size=n,replace=T)
  Data[,"X6"] <- sample(c(0,1),size=n,replace=T)


  true_fail_dists <- function(Dist = "WI"){
    if(Dist == "WI"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(shape = 2, scale = 2.0)
    }else if(Dist == "WD"){
      rdist <- rweibull
      pdist <- pweibull
      params_list <- list(shape = 0.9, scale = 2.0)
    }else if(Dist == "Lgn"){
      rdist <- rlnorm
      pdist <- plnorm
      params_list <- list(meanlog = 1.4, sdlog = 0.4)
    }
    return(list(rdist = rdist,
                pdist = pdist,
                params_list = params_list))
  }

  Dist_info <- true_fail_dists(Dist)

  dat <- LBRC_sampling(n = n,
                       true_fail = Dist_info$rdist,
                       true_params = Dist_info$params_list,
                       ksi = ksi,
                       exp_cens_rate = exp_cens_rate)

  Data$A <- dat[,"left_trunc_time"]
  Data$Y <- dat[,"event_time"]
  Data$event <- dat[,"event"]

  Data$id <- 1:n

  return(Data)
}
