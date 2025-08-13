rm(list=ls())

library(partykit)
library(survival)
library(rsample)

setwd("C:/~")
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

library(Rcpp)
sourceCpp("./methods/vardiCpp.cpp")


Pred_funct <- function(Dist          = "WI",
                       model         = "tree",
                       cov_set_num   = 10,
                       n             = 100,
                       mean_cens     = 20,
                       exp_cens_rate = 1e-100,
                       ksi           = 500,
                       M             = 500)
{
  # All results are stored in the results directory
  base_location <- "./results_tmp"
  # base_location <- "./results"
  final_location <- paste0(base_location,"/",model,"/",Dist,"/N",n)

  if (!dir.exists(final_location)) {
    dir.create(final_location, recursive = TRUE, showWarnings = FALSE)
  }

  setwd(final_location)

  string = paste0("LBRC_DIST_",Dist,
                  "_MODEL_"   ,model,
                  "_P"        ,cov_set_num * 3,
                  "_N"        ,n,
                  "_C"        ,mean_cens,
                  "_sim"      ,M)

  loaded <- tryCatch({
    load(file = string)
    F
  }, error = function(e) {
    T
  })

  if(loaded){
    M_init <- 1

    RES <- list()

    # Tree recovery rate (RR)
    RES$RR$LTRCctree  <- rep(NA,M)
    RES$RR$LBRCctreeC <- rep(NA,M)
    RES$RR$LBRCctreeF <- rep(NA,M)

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

    # learning(tuning) time of forests
    RES$TM_LRN$LTRCcforest   <- rep(NA,M)
    RES$TM_LRN$LBRCcforestCC <- rep(NA,M)
    RES$TM_LRN$LBRCcforestFF <- rep(NA,M)

    # prediction time
    RES$TM_PRD$LTRCcforest   <- rep(NA,M)
    RES$TM_PRD$LBRCcforestCC <- rep(NA,M)
    RES$TM_PRD$LBRCcforestFF <- rep(NA,M)

    # tuned mtry
    RES$opt_mtry$LTRCcforest   <- rep(NA,M)
    RES$opt_mtry$LBRCcforestCC <- rep(NA,M)
    RES$opt_mtry$LBRCcforestFF <- rep(NA,M)

    # mtry that yields minimum L2 difference
    RES$min_mtry$LTRCcforest   <- rep(NA,M)
    RES$min_mtry$LBRCcforestCC <- rep(NA,M)
    RES$min_mtry$LBRCcforestFF <- rep(NA,M)

    # mean censoring rate
    RES$mean_cens_rate <- c()
  }else{
    M_init <- sum(!is.na(RES$L2$LBRCcforestFF$mtryOPT))
    if(M_init>=M){
      cat(string,"Simulation is done!","\n")
      return(RES)
    }
    M_init <- M_init + 1
  }

  # tuning parameter for LBRC-CIT/CIF-F
  vardi_tune <- list(eps=1e-7,max_iter=20)

  # default settings for cforest algorithm
  ntree <- 100
  control = partykit::ctree_control(teststat  = "quad",
                                    testtype  = "Univ",
                                    minsplit  = max(ceiling(sqrt(n)), 20),
                                    minbucket = max(ceiling(sqrt(n)), 7),
                                    saveinfo  = FALSE)

  for(mm in M_init:M){
    ## create the simulation dataset
    if(model %in% c("linear", "nonlinear")){
      nonlinear <- ifelse(model=="nonlinear", TRUE, FALSE)

      # train data
      set.seed(mm)
      DATA <- LBRC.generate_PH_nPH(n             = n,
                                   Dist          = Dist,
                                   exp_cens_rate = exp_cens_rate,
                                   ksi           = ksi,
                                   cov_set_num   = cov_set_num,
                                   nonlinear     = nonlinear)$Data
      # test data
      set.seed(mm+1000)
      TEST_list <- LBRC.generate_PH_nPH(n             = n,
                                        Dist          = Dist,
                                        exp_cens_rate = exp_cens_rate,
                                        ksi           = ksi,
                                        cov_set_num   = cov_set_num,
                                        nonlinear     = nonlinear,
                                        Test_mode     = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv
    }else{
      # train data
      set.seed(mm)
      DATA <- LBRC.generate_tree(n             = n,
                                 Dist          = Dist,
                                 exp_cens_rate = exp_cens_rate,
                                 ksi           = ksi,
                                 cov_set_num   = cov_set_num)$Data
      # test data
      set.seed(mm+1000)
      TEST_list <- LBRC.generate_tree(n             = n,
                                      Dist          = Dist,
                                      exp_cens_rate = exp_cens_rate,
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
    cat("Data mean censoring rate : ",round(mean(RES$mean_cens_rate),2)*100,"%","\n\n",sep="")
    cat("** ", mm, " th L2 distance **","\n",sep="")
    ############################################################################
    ######################### One sample estimator #############################
    ############################################################################
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

    ############################################################################
    ########################### Trees #########################################
    ############################################################################
    if(is.na(RES$L2$LTRCctree[mm])){
      ##################--------- LTRC-CIT --------##################
      obj <- lbrccit(formula = formula, data = DATA,
                     perm_test_est = 'KM')
      if(model == "tree") RES$RR$LTRCctree[mm] <- recoverTree(obj)
      pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                               pred_surv_est = 'KM')
      RES$L2$LTRCctree[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LTRC-CIT       ", "L2 : ", round(RES$L2$LTRCctree[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCctreeCC[mm])){
      ##################--------- LBRC-CIT-(C,C) --------##################
      obj <- lbrccit(formula = formula, data = DATA,
                     perm_test_est = 'MCLE')
      if(model == "tree") RES$RR$LBRCctreeC[mm] <- recoverTree(obj)
      pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                               pred_surv_est = 'MCLE')
      RES$L2$LBRCctreeCC[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LBRC-CIT-C     ", "L2 : ", round(RES$L2$LBRCctreeCC[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCctreeFF[mm])){
      ##################--------- LBRC-CIT-(F,F) --------##################
      obj <- lbrccit(formula = formula, data = DATA,
                     perm_test_est="MFLE", perm_test_args = vardi_tune)
      if(model == "tree") RES$RR$LBRCctreeF[mm] <- recoverTree(obj)
      pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                               pred_surv_est = 'MFLE', pred_surv_args = vardi_tune)
      RES$L2$LBRCctreeFF[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm, " th ", "LBRC-CIT-F     ", "L2 : ", round(RES$L2$LBRCctreeFF[mm],5)*10^5, '\n', sep="")
    }
    ######################### # Activate this code block only for efficiency gain analysis of LBRC-CITs
    # if(is.na(RES$L2$LBRCctreeCF[mm])){
    #   ##################--------- LBRC-CIT-(C,F) --------##################
    #   obj <- lbrccit(formula = formula, data = DATA,
    #                  perm_test_est = 'MCLE')
    #   pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
    #                            pred_surv_est = 'MFLE', pred_surv_args = vardi_tune)
    #   RES$L2$LBRCctreeCF[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    #   cat(mm, " th ", "LBRC-CIT-(C,F) ", "L2 : ", round(RES$L2$LBRCctreeCF[mm],5)*10^5, '\n', sep="")
    # }
    # if(is.na(RES$L2$LBRCctreeFC[mm])){
    #   ##################--------- LBRC-CIT-(F,C) --------##################
    #   obj <- lbrccit(formula = formula, data = DATA,
    #                  perm_test_est="MFLE", perm_test_args = vardi_tune)
    #   pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
    #                            pred_surv_est = 'MCLE')
    #   RES$L2$LBRCctreeFC[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    #   cat(mm, " th ", "LBRC-CIT-(F,C) ", "L2 : ", round(RES$L2$LBRCctreeFC[mm],5)*10^5, '\n', sep="")
    # }

    ###############################################################################
    ################################ cforests #####################################
    ############################################################################
    if(is.na(RES$L2$LTRCcforest$mtryOPT[mm])){
      ##################--------- LTRC-CIF tuned --------##################
      RES$TM_LRN$LTRCcforest[mm] <- system.time({
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=NULL,
                       perm_test_est = "KM",
                       pred_surv_est = "KM",
                       ntree = ntree, control = control, trace = F)
        # Yields identical results to the code below
        # obj <- LTRCforests::ltrccif(formula = formula, data = DATA,
        #                              mtry=length(covariates),
        #                              ntree = ntree, control = control, trace = T)
      })[3]
      RES$TM_PRD$LTRCcforest[mm] <- system.time({
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
      RES$TM_LRN$LBRCcforestCC[mm] <- system.time({
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=NULL,
                       perm_test_est = "MCLE",
                       pred_surv_est = "MCLE",
                       ntree = ntree, control = control, trace = F)
      })[3]
      RES$TM_PRD$LBRCcforestCC[mm] <- system.time({
        pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
                                 pred_surv_est = 'MCLE')
      })[3]
      RES$opt_mtry$LBRCcforestCC[mm] <- obj$mtry
      RES$L2$LBRCcforestCC$mtryOPT[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm," th ", "LBRC-CIF-C     ", "L2 : ", round(RES$L2$LBRCcforestCC$mtryOPT[mm],5)*10^5, '\n', sep="")
    }
    if(is.na(RES$L2$LBRCcforestFF$mtryOPT[mm])){
      ##################--------- LBRC-CIF-F tuned --------##################
      RES$TM_LRN$LBRCcforestFF[mm] <- system.time({
        obj <- lbrccif(formula = formula, data = DATA,
                       mtry=NULL,
                       perm_test_est = "MFLE", perm_test_args = vardi_tune,
                       pred_surv_est = "MFLE", pred_surv_args = vardi_tune,
                       ntree = ntree, control = control, trace = F)
      })[3]
      RES$TM_PRD$LBRCcforestFF[mm] <-

      RES$opt_mtry$LBRCcforestFF[mm] <- obj$mtry
      RES$L2$LBRCcforestFF$mtryOPT[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
      cat(mm," th ", "LBRC-CIF-F     ", "L2 : ", round(RES$L2$LBRCcforestFF$mtryOPT[mm],5)*10^5, '\n', sep="")
    }


    ################## Validating the tuning procedure used for forests #########################
    ################## Activate this code block only for tuning procedure validation
    # mtry_candidates <- c(1,2,3,6,12,24,30)
    # if(is.na(RES[["L2"]][["LTRCcforest"]][["mtryMIN"]][mm])){
    #   cf_L2s <- rep(NA,length(mtry_candidates))
    #   for(i in 1:length(mtry_candidates)){
    #     m <- mtry_candidates[i]
    #     tmp_str <- paste0("mtry",m)
    #     obj <- lbrccif(formula = formula, data = DATA,
    #                    mtry=m,
    #                    perm_test_est = "KM",
    #                    pred_surv_est = "KM",
    #                    ntree = ntree, control = control, trace = T)
    #     pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
    #                              pred_surv_est = 'KM')
    #     RES[["L2"]][["LTRCcforest"]][[tmp_str]][mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    #     cf_L2s[i] <- RES[["L2"]][["LTRCcforest"]][[tmp_str]][mm]
    #     cat("LTRC-CIF",tmp_str,":",round(RES[["L2"]][["LTRCcforest"]][[tmp_str]][mm],5)*10^5, '\n')
    #   }
    #   RES$L2$LTRCcforest$mtryMIN[mm] <- min(cf_L2s)
    #   RES$min_mtry$LTRCcforest[mm] <- mtry_candidates[which.min(cf_L2s)]
    # }
    # if(is.na(RES[["L2"]][["LBRCcforestCC"]][["mtryMIN"]][mm])){
    #   cf_L2s <- rep(NA,length(mtry_candidates))
    #   for(i in 1:length(mtry_candidates)){
    #     m <- mtry_candidates[i]
    #     tmp_str <- paste0("mtry",m)
    #     obj <- lbrccif(formula = formula, data = DATA,
    #                    mtry=m,
    #                    perm_test_est = "MCLE",
    #                    pred_surv_est = "MCLE",
    #                    ntree = ntree, control = control, trace = T)
    #     pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
    #                              pred_surv_est = 'MCLE')
    #     RES[["L2"]][["LBRCcforestCC"]][[tmp_str]][mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    #     cf_L2s[i] <- RES[["L2"]][["LBRCcforestCC"]][[tmp_str]][mm]
    #     cat("LBRC-CIF-C",tmp_str,":",round(RES[["L2"]][["LBRCcforestCC"]][[tmp_str]][mm],5)*10^5, '\n')
    #   }
    #   RES$L2$LBRCcforestCC$mtryMIN[mm] <- min(cf_L2s)
    #   RES$min_mtry$LBRCcforestCC[mm] <- mtry_candidates[which.min(cf_L2s)]
    # }
    # if(is.na(RES[["L2"]][["LBRCcforestFF"]][["mtryMIN"]][mm])){
    #   cf_L2s <- rep(NA,length(mtry_candidates))
    #   for(i in 1:length(mtry_candidates)){
    #     m <- mtry_candidates[i]
    #     tmp_str <- paste0("mtry",m)
    #     obj <- lbrccif(formula = formula, data = DATA,
    #                    mtry=m,
    #                    perm_test_est = "MFLE",
    #                    pred_surv_est = "MFLE",
    #                    ntree = ntree, control = control, trace = T)
    #     pred <- predictProb_LBRC(object = obj, newdata = TEST, time.eval = time.uniq,
    #                              pred_surv_est = 'MFLE')
    #     RES[["L2"]][["LBRCcforestFF"]][[tmp_str]][mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    #     cf_L2s[i] <- RES[["L2"]][["LBRCcforestFF"]][[tmp_str]][mm]
    #     cat("LBRC-CIF-F",tmp_str,":",round(RES[["L2"]][["LBRCcforestFF"]][[tmp_str]][mm],5)*10^5, '\n')
    #   }
    #   RES$L2$LBRCcforestFF$mtryMIN[mm] <- min(cf_L2s)
    #   RES$min_mtry$LBRCcforestFF[mm] <- mtry_candidates[which.min(cf_L2s)]
    # }



    ############################################################################
    ############################# Cox models ###################################
    ############################################################################
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
    ################# print out the simulation summary #########################
    ############################################################################
    cat("** ", mm, "th mean L2 distance **","\n",sep="")
    # tree recovery
    if(model=='tree'){
      cat("LTRC-CIT          RR :", round(mean(RES$RR$LTRCctree[1:mm],na.rm=T),3), '\n')
      cat("LBRC-CIT-(C,)     RR :", round(mean(RES$RR$LBRCctreeC[1:mm],na.rm=T),3), '\n')
      cat("LBRC-CIT-(F,)     RR :", round(mean(RES$RR$LBRCctreeF[1:mm],na.rm=T),3), '\n\n')
    }

    # Integrated L2 difference
    cat("LTRC_1         mean L2 :",round(mean(RES$L2$KM_LT[1:mm],na.rm=T),5)*10^5,'\n')
    cat("MCLE_1         mean L2 :",round(mean(RES$L2$MCLE[1:mm],na.rm=T),5)*10^5,'\n')
    cat("MFLE_1         mean L2 :",round(mean(RES$L2$MFLE[1:mm],na.rm=T),5)*10^5,'\n')

    cat("LTRC-CIT       mean L2 :",round(mean(RES$L2$LTRCctree[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIT-C     mean L2 :",round(mean(RES$L2$LBRCctreeCC[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIT-F     mean L2 :",round(mean(RES$L2$LBRCctreeFF[1:mm],na.rm=T),5)*10^5,'\n')

    cat("LTRC-CIF       mean L2 :",round(mean(RES[["L2"]][["LTRCcforest"]][["mtryOPT"]][1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIF-C     mean L2 :",round(mean(RES[["L2"]][["LBRCcforestCC"]][["mtryOPT"]][1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-CIF-F     mean L2 :",round(mean(RES[["L2"]][["LBRCcforestFF"]][["mtryOPT"]][1:mm],na.rm=T),5)*10^5,'\n')

    cat("LTRC-COX       mean L2 :",round(mean(RES$L2$LTRCcox[1:mm],na.rm=T),5)*10^5,'\n')
    cat("LBRC-COX       mean L2 :",round(mean(RES$L2$LBRCcox[1:mm],na.rm=T),5)*10^5,'\n')

    cat("LTRC-CIF       mean running time :",round(mean(RES$TM_LRN$LTRCcforest + RES$TM_PRD$LTRCcforest, na.rm=T),2),'\n')
    cat("LBRC-CIF-C     mean running time :",round(mean(RES$TM_LRN$LBRCcforestCC + RES$TM_PRD$LBRCcforestCC ,na.rm=T),2),'\n')
    cat("LBRC-CIF-F     mean running time :",round(mean(RES$TM_LRN$LBRCcforestFF + RES$TM_PRD$LBRCcforestFF ,na.rm=T),2),'\n')

    cat("=========================================","\n")
    cat("=========================================","\n\n\n")

    RES.name <- sprintf("LBRC_DIST_%s_MODEL_%s_P%1.0f_N%1.0f_C%1.0f_sim%1.0f",
                        Dist,model,cov_set_num*3,n,mean_cens,M)
    save(RES, file=RES.name)
  }
  return(RES)
}



#################################################################
####################### simulation ##############################
M <- 500            # number of simulation replicates
ksi <- 500          # upper bound for truncation time; set large for the length-biased sampling assumption
cov_set_num <- 10   # number of covariate sets; each set has three types (continuous, ordinal, binary)
ns <- c(100,200,400)  # sample sizes
mean_censs <- c(20,50) # target censoring rates (%)
# exp_cens_rate: rate parameter of the exponential distribution for right-censoring

# =============================== tree ===============================
Dist = "WI" ; model <- "tree"     ; exp_cens_rates <- c(0.084,0.34)
# Dist = "WD" ; model <- "tree"     ; exp_cens_rates <- c(0.08,0.36)
# Dist = "Lgn"; model <- "tree"     ; exp_cens_rates <- c(0.092,0.35)
# Dist = "Bat"; model <- "tree"     ; exp_cens_rates <- c(0.09,0.39)
# =============================== linear =============================
# Dist = "WI" ; model <- "linear"   ; exp_cens_rates <- c(0.14,0.65)
# Dist = "WD" ; model <- "linear"   ; exp_cens_rates <- c(0.12,0.63)
# =============================== nonlinear ==========================
# Dist = "WI" ; model <- "nonlinear"; exp_cens_rates <- c(0.13,0.6)
# Dist = "WD" ; model <- "nonlinear"; exp_cens_rates <- c(0.13,0.63)

for(i in 1:length(ns)){
  n <- ns[i]
  for(c in c(1,2)){
    exp_cens_rate <- exp_cens_rates[c]
    mean_cens     <- mean_censs[c]
    tmp <- Pred_funct(n             = n,
                      Dist          = Dist,
                      model         = model,
                      cov_set_num   = cov_set_num,
                      ksi           = ksi,
                      exp_cens_rate = exp_cens_rate,
                      mean_cens     = mean_cens,
                      M             = M)
  }
}


