library(partykit)
library(survival)
library(rsample)
library(Rcpp)
source("./methods/models_lbrc.R")
source("./methods/plot_lbrc.R")
source("./methods/predictProb_lbrc.R")
source("./methods/sbrier_lbrc.R")
source("./methods/tune.lbrccif.R")
source("./metric/RecoveryRate.R")
source("./metric/L2Distance.R")
source("./data/bathtub.R")
source("./data/generate_LBRC_tree_linear_nonlinear_data.R")
source("./data/LBRC_sampling.R")
sourceCpp("./methods/vardiCpp.cpp")


predict_L2 <- function(model       = "tree",
                       Dist        = "WI",
                       cov_set_num = 10,
                       n           = 100,
                       cens_rate   = 20,
                       ksi         = 500,
                       M           = 500,
                       working_dir = NULL)
{
  # All results are stored in the results directory
  result_location <- paste0(working_dir,"/results_intermediate","/",model,"/",Dist,"/N",n)

  if (!dir.exists(result_location)) {
    dir.create(result_location, recursive = TRUE, showWarnings = FALSE)
  }

  setwd(result_location)

  string = paste0("LBRC_DIST_",Dist,
                  "_MODEL_"   ,model,
                  "_P"        ,cov_set_num * 3,
                  "_N"        ,n,
                  "_C"        ,cens_rate)

  file_missing <- tryCatch({
    load(file = string)
    F
  }, error = function(e) {
    T
  })

  if(file_missing){
    RES <- list()

    # Tree recovery rate (RR)
    RES$RR$LTRCctree   <- rep(NA,M)
    RES$RR$LBRCctreeC  <- rep(NA,M)
    RES$RR$LBRCctreeF  <- rep(NA,M)

    # Integrated L2 error
    RES$L2$KM_LT       <- rep(NA,M)
    RES$L2$MCLE        <- rep(NA,M)
    RES$L2$MFLE        <- rep(NA,M)
    RES$L2$LTRCctree   <- rep(NA,M)
    RES$L2$LBRCctreeCC <- rep(NA,M)
    RES$L2$LBRCctreeFF <- rep(NA,M)
    RES$L2$LBRCctreeCF <- rep(NA,M)
    RES$L2$LBRCctreeFC <- rep(NA,M)
    mtry_candidates <- c(1,2,3,6,12,24,30)
    for(m in mtry_candidates){
      tmp_str <- paste0("mtry",m)
      RES[["L2"]][["LTRCcforest"]][[tmp_str]]   <- rep(NA,M)
      RES[["L2"]][["LBRCcforestCC"]][[tmp_str]] <- rep(NA,M)
      RES[["L2"]][["LBRCcforestFF"]][[tmp_str]] <- rep(NA,M)
    }
    RES[["L2"]][["LTRCcforest"]][["mtryOPT"]]   <- rep(NA,M)
    RES[["L2"]][["LBRCcforestCC"]][["mtryOPT"]] <- rep(NA,M)
    RES[["L2"]][["LBRCcforestFF"]][["mtryOPT"]] <- rep(NA,M)
    RES[["L2"]][["LTRCcforest"]][["mtryMIN"]]   <- rep(NA,M)
    RES[["L2"]][["LBRCcforestCC"]][["mtryMIN"]] <- rep(NA,M)
    RES[["L2"]][["LBRCcforestFF"]][["mtryMIN"]] <- rep(NA,M)
    RES$L2$LTRCcox     <- rep(NA,M)
    RES$L2$LBRCcox     <- rep(NA,M)

    # Running time (RN_TM) of trees and forests
    RES$RN_TM$LTRCctree     <- rep(NA,M)
    RES$RN_TM$LBRCctreeCC   <- rep(NA,M)
    RES$RN_TM$LBRCctreeFF   <- rep(NA,M)
    RES$RN_TM$LTRCcforest   <- rep(NA,M)
    RES$RN_TM$LBRCcforestCC <- rep(NA,M)
    RES$RN_TM$LBRCcforestFF <- rep(NA,M)

    # Tuned mtry
    RES$opt_mtry$LTRCcforest   <- rep(NA,M)
    RES$opt_mtry$LBRCcforestCC <- rep(NA,M)
    RES$opt_mtry$LBRCcforestFF <- rep(NA,M)

    # mtry that yields minimum L2 difference
    RES$min_mtry$LTRCcforest   <- rep(NA,M)
    RES$min_mtry$LBRCcforestCC <- rep(NA,M)
    RES$min_mtry$LBRCcforestFF <- rep(NA,M)

    # Mean censoring rate for current simulation setting
    RES$mean_cens_rate <- c()
  }

  M_init <- sum(!is.na(RES$L2$LBRCcox))
  if(M_init >= M){
    cat(string,"Simulation is done for",M,"!","\n")
    return(RES)
  }
  M_init <- M_init + 1

  # Tuning parameter for LBRC-CIT/CIF-F
  vardi_tune <- list(eps=1e-7,max_iter=20)

  # Default settings for cforest algorithm
  ntree <- 100
  control = partykit::ctree_control(teststat  = "quad",
                                    testtype  = "Univ",
                                    minsplit  = max(ceiling(sqrt(n)), 20),
                                    minbucket = max(ceiling(sqrt(n)), 7),
                                    saveinfo  = FALSE)
  for(mm in M_init:M){
    ## Create the simulation data
    if(model %in% c("linear", "nonlinear")){ # generate (non)linear data
      nonlinear <- ifelse(model=="nonlinear", TRUE, FALSE)

      # Generate train data
      set.seed(mm)
      DATA <- LBRC.generate_PH_nPH(n             = n,
                                   Dist          = Dist,
                                   cens_rate     = cens_rate,
                                   ksi           = ksi,
                                   cov_set_num   = cov_set_num,
                                   nonlinear     = nonlinear)$Data
      # Generate test data
      set.seed(mm + 1000)
      TEST_list <- LBRC.generate_PH_nPH(n             = n,
                                        Dist          = Dist,
                                        cens_rate     = cens_rate,
                                        ksi           = ksi,
                                        cov_set_num   = cov_set_num,
                                        nonlinear     = nonlinear,
                                        Test_mode     = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv
    }else{ # generate tree structured data
      # Generate train data
      set.seed(mm)
      DATA <- LBRC.generate_tree(n             = n,
                                 Dist          = Dist,
                                 cens_rate     = cens_rate,
                                 ksi           = ksi,
                                 cov_set_num   = cov_set_num)$Data
      # Generate test data
      set.seed(mm + 1000)
      TEST_list <- LBRC.generate_tree(n             = n,
                                      Dist          = Dist,
                                      cens_rate     = cens_rate,
                                      ksi           = ksi,
                                      cov_set_num   = cov_set_num,
                                      Test_mode     = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv
    }

    covariates <- setdiff(colnames(DATA), c("A","Y","event","id"))
    formula <- as.formula(paste("Surv(A, Y, event) ~", paste(covariates, collapse = "+")))

    RES$mean_cens_rate[mm] <- mean(1-DATA$event)

    cat("=========================================","\n")
    cat("=========================================","\n")
    cat("** ", mm, " th data successfully generated ! **","\n",sep="")
    cat("Data Structure           : ",model,"\n",sep="")
    cat("Data Distribution        : ",Dist,"\n",sep="")
    cat("Data Sample size         : ",n,"\n",sep="")
    cat("Data mean censoring rate : ",round(mean(RES$mean_cens_rate),2)*100,"%","\n",sep="")
    cat("** ", mm, " th L2 distance **","\n",sep="")

    ##########################################################################################
    ##########################################################################################
    ##########                                                                      ##########
    ##########                        One sample estimator                          ##########
    ##########                                                                      ##########
    ##########################################################################################
    ##########################################################################################
    y_mat <- as.matrix(DATA[,c("A","Y","event")])
    if(is.na(RES$L2$KM_LT[mm])){
      ##################--------- LTRC stump --------##################
      obj <- .pred_Surv_nolog(y=Surv(y_mat[,1],y_mat[,2],y_mat[,3]), w=rep(1,dim(DATA)[1]))
      pred <- lapply(1:n, function(i) obj)
      pred <- predict_onesample_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
      RES$L2$KM_LT[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LTRC-1         ", "L2 : ", round(RES$L2$KM_LT[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$MCLE[mm])){
      ##################--------- MCLE stump --------##################
      obj <- .pred_Surv_LBRC_MCLE(y=y_mat, w=rep(1,dim(DATA)[1]))
      pred <- lapply(1:n, function(i) obj)
      pred <- predict_onesample_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
      RES$L2$MCLE[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "MCLE-1         ", "L2 : ", round(RES$L2$MCLE[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$MFLE[mm])){
      ##################--------- MFLE stump --------##################
      obj <- .pred_Surv_LBRC_MFLE(y=y_mat, w=rep(1,dim(DATA)[1]),
                                  eps=vardi_tune$eps, max_iter=vardi_tune$max_iter)
      pred <- lapply(1:n, function(i) obj)
      pred <- predict_onesample_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
      RES$L2$MFLE[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "MFLE-1         ", "L2 : ", round(RES$L2$MFLE[mm],5)*10^5, '\n', sep="")
    }

    ##########################################################################################
    ##########################################################################################
    ##########                                                                      ##########
    ##########                      Conditional inference tree                      ##########
    ##########                                                                      ##########
    ##########################################################################################
    ##########################################################################################
    if(is.na(RES$L2$LTRCctree[mm])){
      ##################--------- LTRC-CIT --------##################
      RES$RN_TM$LTRCctree[mm] <- system.time({
        obj <- lbrccit(formula = formula, data = DATA,
                       perm_test_est = 'KM')
        # Yields identical results to the code below
        # obj <- LTRCtrees::LTRCIT(Formula = formula, data = DATA)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'KM')
      })[3]
      if(model == "tree") RES$RR$LTRCctree[mm] <- recoverTree(obj)
      RES$L2$LTRCctree[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LTRC-CIT       ", "L2 : ", round(RES$L2$LTRCctree[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCctreeCC[mm])){
      ##################--------- LBRC-CIT-(C,C) --------##################
      RES$RN_TM$LBRCctreeCC[mm] <- system.time({
        obj <- lbrccit(formula = formula, data = DATA,
                       perm_test_est = 'MCLE')
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MCLE')
      })[3]
      if(model == "tree") RES$RR$LBRCctreeC[mm] <- recoverTree(obj)
      RES$L2$LBRCctreeCC[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LBRC-CIT-(C,C) ", "L2 : ", round(RES$L2$LBRCctreeCC[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCctreeFF[mm])){
      ##################--------- LBRC-CIT-(F,F) --------##################
      RES$RN_TM$LBRCctreeFF[mm] <- system.time({
        obj <- lbrccit(formula = formula, data = DATA,
                       perm_test_est="MFLE", perm_test_args = vardi_tune)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MFLE', pred_surv_args = vardi_tune)
      })[3]
      if(model == "tree") RES$RR$LBRCctreeF[mm] <- recoverTree(obj)
      RES$L2$LBRCctreeFF[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LBRC-CIT-(F,F) ", "L2 : ", round(RES$L2$LBRCctreeFF[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCctreeCF[mm])){
      ##################--------- LBRC-CIT-(C,F) --------##################
      obj <- lbrccit(formula = formula, data = DATA,
                     perm_test_est = 'MCLE')
      pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                               pred_surv_est = 'MFLE', pred_surv_args = vardi_tune)
      RES$L2$LBRCctreeCF[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LBRC-CIT-(C,F) ", "L2 : ", round(RES$L2$LBRCctreeCF[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCctreeFC[mm])){
      ##################--------- LBRC-CIT-(F,C) --------##################
      obj <- lbrccit(formula = formula, data = DATA,
                     perm_test_est="MFLE", perm_test_args = vardi_tune)
      pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                               pred_surv_est = 'MCLE')
      RES$L2$LBRCctreeFC[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LBRC-CIT-(F,C) ", "L2 : ", round(RES$L2$LBRCctreeFC[mm],5)*10^5, '\n', sep="")
    }

    ##########################################################################################
    ##########################################################################################
    ##########                                                                      ##########
    ##########                      Conditional inference forest                    ##########
    ##########                                                                      ##########
    ##########################################################################################
    ##########################################################################################
    if(is.na(RES$L2$LTRCcforest$mtryOPT[mm])){
      ##################--------- LTRC-CIF tuned --------##################
      RES$RN_TM$LTRCcforest[mm] <- system.time({
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=NULL,
                       perm_test_est = "KM",
                       pred_surv_est = "KM",
                       ntree = ntree, control = control, trace = F)
        # Yields identical results to the code below
        # obj <- LTRCforests::ltrccif(formula = formula, data = DATA,
        #                              mtry=length(covariates),
        #                              ntree = ntree, control = control, trace = T)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'KM')
        # Yields identical results to the code below
        # pred <- LTRCforests::predictProb(object = obj, newdata = TEST, time.eval = time.uniq)
      })[3]
      RES$opt_mtry$LTRCcforest[mm] <- obj$mtry
      RES$L2$LTRCcforest$mtryOPT[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LTRC-CIF       ", "L2 : ", round(RES$L2$LTRCcforest$mtryOPT[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCcforestCC$mtryOPT[mm])){
      ##################--------- LBRC-CIF-C tuned --------##################
      RES$RN_TM$LBRCcforestCC[mm] <- system.time({
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=NULL,
                       perm_test_est = "MCLE",
                       pred_surv_est = "MCLE",
                       ntree = ntree, control = control, trace = F)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MCLE')
      })[3]
      RES$opt_mtry$LBRCcforestCC[mm] <- obj$mtry
      RES$L2$LBRCcforestCC$mtryOPT[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm," th ", "LBRC-CIF-C     ", "L2 : ", round(RES$L2$LBRCcforestCC$mtryOPT[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCcforestFF$mtryOPT[mm])){
      ##################--------- LBRC-CIF-F tuned --------##################
      RES$RN_TM$LBRCcforestFF[mm] <- system.time({
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=NULL,
                       perm_test_est = "MFLE", perm_test_args = vardi_tune,
                       pred_surv_est = "MFLE", pred_surv_args = vardi_tune,
                       ntree = ntree, control = control, trace = F)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MFLE')
      })[3]
      RES$opt_mtry$LBRCcforestFF[mm] <- obj$mtry
      RES$L2$LBRCcforestFF$mtryOPT[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm," th ", "LBRC-CIF-F     ", "L2 : ", round(RES$L2$LBRCcforestFF$mtryOPT[mm],5)*10^5, '\n', sep="")
    }

    ##########################################################################################
    ##########################################################################################
    ##########                                                                      ##########
    ##########                               Cox model                              ##########
    ##########                                                                      ##########
    ##########################################################################################
    ##########################################################################################
    if(is.na(RES$L2$LTRCcox[mm])){
      ##################--------- L2 for LTRC-COX --------##################
      obj <- coxph(formula = formula, data = DATA)
      pred <- survfit(obj, newdata = TEST)
      pred <- lapply(1:n, function(i) pred[i])
      pred <- predict_onesample_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
      RES$L2$LTRCcox[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm," th ", "LTRC-COX       ", "L2 : ", round(RES$L2$LTRCcox[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCcox[mm])){
      ##################--------- L2 for LBRC-COX --------##################
      # from supplement R-code https://pmc.ncbi.nlm.nih.gov/articles/PMC3758493/#SM
      CPDATA <- data.frame( cbind(
        A = c(DATA$A, (DATA$Y-DATA$A)[DATA$event == 1]),
        Y = c(DATA$Y,  DATA$Y[DATA$event == 1]),
        event = c(DATA$event,  DATA$event[DATA$event == 1])
      ))
      for(cov in covariates){
        CPDATA[,cov] <- c(DATA[,cov], DATA[DATA$event==1,cov])
      }
      cpfit <- coxph(formula, data = CPDATA)
      pred <- survfit(cpfit, newdata = TEST)
      pred <- lapply(1:n, function(i) pred[i])
      pred <- predict_onesample_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
      RES$L2$LBRCcox[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm," th ", "LBRC-COX       ", "L2 : ", round(RES$L2$LBRCcox[mm],5)*10^5, '\n', sep="")
    }

    ############################################################################
    ##################### print out the simulation summary #####################
    ############################################################################
    cat("\n", "** ", mm, "th mean results **","\n",sep="")
    # Mean - tree recovery rate
    if(model=='tree'){
      cat("LTRC-CIT          RR :", round(mean(RES$RR$LTRCctree[1:mm],na.rm=T),3),'\n')
      cat("LBRC-CIT-(C,)     RR :", round(mean(RES$RR$LBRCctreeC[1:mm],na.rm=T),3),'\n')
      cat("LBRC-CIT-(F,)     RR :", round(mean(RES$RR$LBRCctreeF[1:mm],na.rm=T),3),'\n')
    }

    # Mean - integrated L2 difference
    cat("LTRC-stump     mean L2 :",round(mean(RES$L2$KM_LT[1:mm],na.rm=T),5)*10^5,'\n')
    cat("MCLE-stump     mean L2 :",round(mean(RES$L2$MCLE[1:mm],na.rm=T),5)*10^5,'\n')
    cat("MFLE-stump     mean L2 :",round(mean(RES$L2$MFLE[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LTRC-CIT       mean L2 :",round(mean(RES$L2$LTRCctree[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIT-C     mean L2 :",round(mean(RES$L2$LBRCctreeCC[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIT-F     mean L2 :",round(mean(RES$L2$LBRCctreeFF[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LTRC-CIF       mean L2 :",round(mean(RES[["L2"]][["LTRCcforest"]][["mtryOPT"]][1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIF-C     mean L2 :",round(mean(RES[["L2"]][["LBRCcforestCC"]][["mtryOPT"]][1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIF-F     mean L2 :",round(mean(RES[["L2"]][["LBRCcforestFF"]][["mtryOPT"]][1:mm],na.rm=T),5)*10^5,'\n')
    cat("LTRC-COX       mean L2 :",round(mean(RES$L2$LTRCcox[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-COX       mean L2 :",round(mean(RES$L2$LBRCcox[1:mm],na.rm=T),5)*10^5,'\n')
    cat("=========================================","\n")
    cat("=========================================","\n\n\n")

    RES.name <- sprintf("LBRC_DIST_%s_MODEL_%s_P%1.0f_N%1.0f_C%1.0f",
                        Dist,model,cov_set_num*3,n,cens_rate)
    save(RES, file=RES.name)
  }
}

validate_OOB_tuning <- function(model       = "tree",
                                Dist        = "WI",
                                cov_set_num = 10,
                                n           = 100,
                                cens_rate   = 20,
                                ksi         = 500,
                                M           = 500,
                                working_dir = NULL)
{
  # All results are stored in the results directory
  result_location <- paste0(working_dir,"/results_intermediate","/",model,"/",Dist,"/N",n)

  if (!dir.exists(result_location)) {
    dir.create(result_location, recursive = TRUE, showWarnings = FALSE)
  }

  setwd(result_location)

  string = paste0("LBRC_DIST_",Dist,
                  "_MODEL_"   ,model,
                  "_P"        ,cov_set_num * 3,
                  "_N"        ,n,
                  "_C"        ,cens_rate)

  file_missing <- tryCatch({
    load(file = string)
    F
  }, error = function(e) {
    T
  })

  if(file_missing){
    cat("There is no stored file: run other simulation first!")
    break
  }

  M_init <- sum(!is.na(RES$L2$LBRCcforestFF$mtryMIN))
  if(M_init >= M){
    cat(string,"Simulation is done for",M,"!","\n")
    return(RES)
  }
  M_init <- M_init + 1

  # Tuning parameter for LBRC-CIT/CIF-F
  vardi_tune <- list(eps=1e-7,max_iter=20)

  # Default settings for cforest algorithm
  ntree <- 100
  control = partykit::ctree_control(teststat  = "quad",
                                    testtype  = "Univ",
                                    minsplit  = max(ceiling(sqrt(n)), 20),
                                    minbucket = max(ceiling(sqrt(n)), 7),
                                    saveinfo  = FALSE)

  for(mm in M_init:M){
    ## Create the simulation data
    if(model %in% c("linear", "nonlinear")){ # generate (non)linear data
      nonlinear <- ifelse(model=="nonlinear", TRUE, FALSE)

      # Generate train data
      set.seed(mm)
      DATA <- LBRC.generate_PH_nPH(n             = n,
                                   Dist          = Dist,
                                   cens_rate     = cens_rate,
                                   ksi           = ksi,
                                   cov_set_num   = cov_set_num,
                                   nonlinear     = nonlinear)$Data
      # Generate test data
      set.seed(mm + 1000)
      TEST_list <- LBRC.generate_PH_nPH(n             = n,
                                        Dist          = Dist,
                                        cens_rate     = cens_rate,
                                        ksi           = ksi,
                                        cov_set_num   = cov_set_num,
                                        nonlinear     = nonlinear,
                                        Test_mode     = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv
    }else{ # generate tree structured data
      # Generate train data
      set.seed(mm)
      DATA <- LBRC.generate_tree(n             = n,
                                 Dist          = Dist,
                                 cens_rate     = cens_rate,
                                 ksi           = ksi,
                                 cov_set_num   = cov_set_num)$Data
      # Generate test data
      set.seed(mm + 1000)
      TEST_list <- LBRC.generate_tree(n             = n,
                                      Dist          = Dist,
                                      cens_rate     = cens_rate,
                                      ksi           = ksi,
                                      cov_set_num   = cov_set_num,
                                      Test_mode     = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv
    }

    covariates <- setdiff(colnames(DATA), c("A","Y","event","id"))
    formula <- as.formula(paste("Surv(A, Y, event) ~", paste(covariates, collapse = "+")))

    RES$mean_cens_rate[mm] <- mean(1-DATA$event)

    cat("=========================================","\n")
    cat("=========================================","\n")
    cat("** ", mm, " th data successfully generated ! **","\n",sep="")
    cat("Data Structure           : ",model,"\n",sep="")
    cat("Data Distribution        : ",Dist,"\n",sep="")
    cat("Data Sample size         : ",n,"\n",sep="")
    cat("Data mean censoring rate : ",round(mean(RES$mean_cens_rate),2)*100,"%","\n",sep="")
    cat("** ", mm, " th L2 distance **","\n",sep="")

    ##########################################################################################
    ##########################################################################################
    ##########                                                                      ##########
    ##########          Validation of the tuning procedure used for forests         ##########
    ##########                                                                      ##########
    ##########################################################################################
    ##########################################################################################
    mtry_candidates <- c(1,2,3,6,12,24,30)
    if(is.na(RES[["L2"]][["LTRCcforest"]][["mtryMIN"]][mm])){
      cf_L2s <- rep(NA,length(mtry_candidates))
      for(i in 1:length(mtry_candidates)){
        m <- mtry_candidates[i]
        tmp_str <- paste0("mtry",m)
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=m,
                       perm_test_est = "KM",
                       pred_surv_est = "KM",
                       ntree = ntree, control = control, trace = T)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'KM')
        RES[["L2"]][["LTRCcforest"]][[tmp_str]][mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
        cf_L2s[i] <- RES[["L2"]][["LTRCcforest"]][[tmp_str]][mm]
        cat("LTRC-CIF",tmp_str,":",round(RES[["L2"]][["LTRCcforest"]][[tmp_str]][mm],5)*10^5, '\n')
      }
      RES$L2$LTRCcforest$mtryMIN[mm] <- min(cf_L2s)
      RES$min_mtry$LTRCcforest[mm] <- mtry_candidates[which.min(cf_L2s)]
    }
    if(is.na(RES[["L2"]][["LBRCcforestCC"]][["mtryMIN"]][mm])){
      cf_L2s <- rep(NA,length(mtry_candidates))
      for(i in 1:length(mtry_candidates)){
        m <- mtry_candidates[i]
        tmp_str <- paste0("mtry",m)
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=m,
                       perm_test_est = "MCLE",
                       pred_surv_est = "MCLE",
                       ntree = ntree, control = control, trace = T)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MCLE')
        RES[["L2"]][["LBRCcforestCC"]][[tmp_str]][mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
        cf_L2s[i] <- RES[["L2"]][["LBRCcforestCC"]][[tmp_str]][mm]
        cat("LBRC-CIF-C",tmp_str,":",round(RES[["L2"]][["LBRCcforestCC"]][[tmp_str]][mm],5)*10^5, '\n')
      }
      RES$L2$LBRCcforestCC$mtryMIN[mm] <- min(cf_L2s)
      RES$min_mtry$LBRCcforestCC[mm] <- mtry_candidates[which.min(cf_L2s)]
    }
    if(is.na(RES[["L2"]][["LBRCcforestFF"]][["mtryMIN"]][mm])){
      cf_L2s <- rep(NA,length(mtry_candidates))
      for(i in 1:length(mtry_candidates)){
        m <- mtry_candidates[i]
        tmp_str <- paste0("mtry",m)
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=m,
                       perm_test_est = "MFLE",
                       pred_surv_est = "MFLE",
                       ntree = ntree, control = control, trace = T)
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MFLE')
        RES[["L2"]][["LBRCcforestFF"]][[tmp_str]][mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
        cf_L2s[i] <- RES[["L2"]][["LBRCcforestFF"]][[tmp_str]][mm]
        cat("LBRC-CIF-F",tmp_str,":",round(RES[["L2"]][["LBRCcforestFF"]][[tmp_str]][mm],5)*10^5, '\n')
      }
      RES$L2$LBRCcforestFF$mtryMIN[mm] <- min(cf_L2s)
      RES$min_mtry$LBRCcforestFF[mm] <- mtry_candidates[which.min(cf_L2s)]
    }

    ############################################################################
    ##################### print out the simulation summary #####################
    ############################################################################
    cat("\n", "** ", mm, "th mean results **","\n",sep="")
    cat("=========================================","\n")
    cat("=========================================","\n\n\n")

    RES.name <- sprintf("LBRC_DIST_%s_MODEL_%s_P%1.0f_N%1.0f_C%1.0f",
                        Dist,model,cov_set_num*3,n,cens_rate)
    save(RES, file=RES.name)
  }
}

test_unbiasedness <- function(Dist        = "WI",
                              n           = 200,
                              cens_rate   = 20,
                              ksi         = 500,
                              M           = 10000,
                              working_dir = NULL)
{
  # All results are stored in the results directory
  result_location <- paste0(working_dir,"/results_intermediate/test_unbiasedness")

  if (!dir.exists(result_location)) {
    dir.create(result_location, recursive = TRUE, showWarnings = FALSE)
  }

  setwd(result_location)

  string = paste0("LBRC_UNBIASED_TEST_DIST_",Dist,
                  "_C"                      ,cens_rate)

  file_missing <- tryCatch({
    load(file = string)
    F
  }, error = function(e) {
    T
  })

  if(file_missing){
    RES <- list()

    # Variable chosen at tree root
    RES$RV$LBRCctreeC <- rep(NA,M)
    RES$RV$LBRCctreeF <- rep(NA,M)

    # Mean censoring rate for current simulation setting
    RES$mean_cens_rate <- c()
  }

  M_init <- sum(!is.na(RES$RV$LBRCctreeF))
  if(M_init >= M){
    ctable <- table(RES$RV$LBRCctreeC)
    ftable <- table(RES$RV$LBRCctreeF)
    cp <- round(chisq.test(ctable)$p.value, 4)
    fp <- round(chisq.test(ftable)$p.value, 4)

    cat("=======================================================================","\n")
    cat("=======================================================================","\n")
    cat("** ", M, "th Table of selected variables under mean censoring rate : ", round(mean(RES$mean_cens_rate),2)*100, "% **\n",sep="")
    cat(sprintf("%-16s %5s %5s %5s %5s %5s %5s %-7s\n",
                "", "X1", "X2", "X3", "X4", "X5", "X6", "p-value"))
    cat(sprintf("%-16s %5d %5d %5d %5d %5d %5d %7.4f\n",
                "LBRC-CIT-(C,) :", ctable[1], ctable[2], ctable[3], ctable[4], ctable[5], ctable[6], cp))
    cat(sprintf("%-16s %5d %5d %5d %5d %5d %5d %7.4f\n",
                "LBRC-CIT-(F,) :", ftable[1], ftable[2], ftable[3], ftable[4], ftable[5], ftable[6], fp))
    cat("=======================================================================","\n")
    cat("=======================================================================","\n\n\n")
    return(RES)
  }
  M_init <- M_init + 1

  # important!! This forces splitting
  control = partykit::ctree_control(
    mincriterion = 0,
    stump        = T,
    minsplit     = 2,
    minbucket    = 1
  )

  # Tuning parameter for LBRC-CIT-F
  vardi_tune <- list(eps=1e-7,max_iter=20)

  cat("============================================================================","\n")
  cat("===== Test of unbiasedness under ", Dist, " distribution and ", cens_rate, "%", " rate started! =====", "\n", sep="")
  cat("============================================================================","\n")

  for(mm in M_init:M){
    ## Create the simulation data
    set.seed(mm)
    DATA <- LBRC.generate_one_sample(n         = n,
                                     Dist      = Dist,
                                     cens_rate = cens_rate,
                                     ksi       = ksi)

    covariates <- setdiff(colnames(DATA), c("A","Y","event","id"))
    formula <- as.formula(paste("Surv(A, Y, event) ~", paste(covariates, collapse = "+")))

    RES$mean_cens_rate[mm] <- mean(1-DATA$event)

    ##################--------- LBRC-CIT-(C,) --------##################
    obj <- lbrccit(formula = formula, data = DATA,
                   perm_test_est = 'MCLE',
                   control = control)
    RES$RV$LBRCctreeC[mm] <- root_node_return(obj)

    ##################--------- LBRC-CIT-(F,) --------##################
    obj <- lbrccit(formula = formula, data = DATA,
                   perm_test_est="MFLE", perm_test_args = vardi_tune,
                   control = control)
    RES$RV$LBRCctreeF[mm] <- root_node_return(obj)

    if(mm%%100 == 0 | mm == M){
      ctable <- table(RES$RV$LBRCctreeC)
      ftable <- table(RES$RV$LBRCctreeF)
      cp <- round(chisq.test(ctable)$p.value, 4)
      fp <- round(chisq.test(ftable)$p.value, 4)

      cat("=======================================================================","\n")
      cat("=======================================================================","\n")
      cat("** ", mm, "th Table of selected variables under mean censoring rate : ", round(mean(RES$mean_cens_rate),2)*100, "% **\n",sep="")
      cat(sprintf("%-16s %5s %5s %5s %5s %5s %5s %-7s\n",
                  "", "X1", "X2", "X3", "X4", "X5", "X6", "p-value"))
      cat(sprintf("%-16s %5d %5d %5d %5d %5d %5d %7.4f\n",
                  "LBRC-CIT-(C,) :", ctable[1], ctable[2], ctable[3], ctable[4], ctable[5], ctable[6], cp))
      cat(sprintf("%-16s %5d %5d %5d %5d %5d %5d %7.4f\n",
                  "LBRC-CIT-(F,) :", ftable[1], ftable[2], ftable[3], ftable[4], ftable[5], ftable[6], fp))
      cat("=======================================================================","\n")
      cat("=======================================================================","\n\n\n")
    }

    RES.name <- sprintf("LBRC_UNBIASED_TEST_DIST_%s_C%1.0f",
                        Dist,cens_rate)
    save(RES, file=RES.name)
  }
}

simulate_LBRC_tree_methods <- function(mode = "model comparison",
                                       sim_set_list, M, working_dir){
  for(i in seq_len(length(sim_set_list$ns))){
    for(c in seq_len(length(sim_set_list$cens_rates))){
      if(mode == "model prediction"){
        predict_L2(model       = sim_set_list$model,
                   Dist        = sim_set_list$Dist,
                   cov_set_num = sim_set_list$cov_set_num,
                   n           = sim_set_list$ns[i],
                   cens_rate   = sim_set_list$cens_rates[c],
                   ksi         = sim_set_list$ksi,
                   M           = M,
                   working_dir = working_dir)
      }else if(mode == "OOB tuning validation"){
        validate_OOB_tuning(model       = sim_set_list$model,
                            Dist        = sim_set_list$Dist,
                            cov_set_num = sim_set_list$cov_set_num,
                            n           = sim_set_list$ns[i],
                            cens_rate   = sim_set_list$cens_rates[c],
                            ksi         = sim_set_list$ksi,
                            M           = M,
                            working_dir = working_dir)
      }else if(mode == "test unbiasedness"){
        test_unbiasedness(Dist        = sim_set_list$Dist,
                          n           = sim_set_list$ns[i],
                          cens_rate   = sim_set_list$cens_rates[c],
                          ksi         = sim_set_list$ksi,
                          M           = M,
                          working_dir = working_dir)
      }else{
        cat("Please select the mode first!")
        break
      }
    }
  }
}







