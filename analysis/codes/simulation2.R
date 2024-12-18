rm(list=ls())

setwd("C:/Users/jinwu/Desktop/LBRCtree_forests/pkg/LBRCtreeforests")
devtools::load_all()

# install.packages("C:/Users/jinwu/Desktop/LBRCtree_forests/analysis/utils/CoxPhLb_1.2.0.tar.gz")
library(CoxPhLb)
library(partykit)
library(survival)

setwd("C:/Users/jinwu/Desktop/LBRCtree_forests/analysis")
source("./utils/f1_macro.R")
source("./utils/LossFunct.R")
source("./data/Length_Biased_generation.R")

Pred_funct <- function(n=100, Dist="Exp", model="tree",
                       cov_set_num=2, exp_cens_rate=1e-8, ksi=100, M=100,
                       mean_cens = 20){
  RES.RR = NULL
  RES.RR$RCctree <- rep(NA,M)
  RES.RR$LTRCctree <- rep(NA,M)
  RES.RR$LBRCctree <- rep(NA,M)

  RES.F1 = NULL
  RES.F1$RCctree <- rep(NA,M)
  RES.F1$LTRCctree <- rep(NA,M)
  RES.F1$LBRCctree <- rep(NA,M)

  RES.L2 = NULL
  RES.L2$RCctree <- rep(NA,M)
  RES.L2$LTRCctree <- rep(NA,M)
  RES.L2$LBRCctree <- rep(NA,M)
  RES.L2$RCcforest.D <- rep(NA,M); RES.L2$RCcforest.T <- rep(NA,M)
  RES.L2$LTRCcforest.D <- rep(NA,M); RES.L2$LTRCcforest.T <- rep(NA,M)
  RES.L2$LBRCcforest.D <- rep(NA,M); RES.L2$LBRCcforest.T <- rep(NA,M)
  RES.L2$RCcox <- rep(NA,M)
  RES.L2$LTRCcox <- rep(NA,M)
  RES.L2$LBRCcox.I <- rep(NA,M)
  RES.L2$LBRCcox.II <- rep(NA,M)

  # for recording mean censoring rate over M simulations
  cens_rate_list <- c()

  # default settings for cforest algorithm
  ntree <- 100
  control = partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                    # minsplit = max(ceiling(sqrt(nrow(data))), 20),
                                    # minbucket = max(ceiling(sqrt(nrow(data))), 7),
                                    minsplit = n * 0.15,
                                    minbucket = n * 0.06,
                                    minprob = 0.01,
                                    # maxdepth = 2,
                                    # mincriterion = 0,
                                    saveinfo = FALSE)

  for(mm in 1:M){
    ## create the simulation dataset
    if(model %in% c("linear", "nonlinear")){
      nonlinear <- ifelse(model=="nonlinear", TRUE, FALSE)

      set.seed(mm)
      DATA <- LBRC.generate_PH_nPH(n = n,
                                   Dist = Dist,
                                   exp_cens_rate = exp_cens_rate,
                                   ksi = ksi,
                                   cov_set_num = cov_set_num,
                                   nonlinear = nonlinear)$Data
      set.seed(mm+M)
      TEST_list <- LBRC.generate_PH_nPH(n = n,
                                        Dist = Dist,
                                        exp_cens_rate = exp_cens_rate,
                                        ksi = ksi,
                                        cov_set_num = cov_set_num,
                                        nonlinear = nonlinear,
                                        Test_mode = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv_T_A
    }else{
      set.seed(mm)
      DATA <- LBRC.generate_tree(n = n,
                                 Dist = Dist,
                                 exp_cens_rate = exp_cens_rate,
                                 ksi = ksi,
                                 cov_set_num = cov_set_num)$Data
      set.seed(mm+M)
      TEST_list <- LBRC.generate_tree(n = n,
                                      Dist = Dist,
                                      exp_cens_rate = exp_cens_rate,
                                      ksi = ksi,
                                      cov_set_num = cov_set_num,
                                      Test_mode = TRUE)
      TEST <- TEST_list$Data
      time.uniq <- TEST_list$time.uniq
      TEST.true <- TEST_list$true_surv_T_A
      # essential labels to calculate RR and F1
      true_vars <- c("X1","X2",'X3')
      tr_trmnd <- TEST_list$tr_trmnd
    }

    covariates <- setdiff(colnames(DATA), c("A","Y","event","id"))
    formula <- as.formula(paste("Surv(A, Y, event) ~", paste(covariates, collapse = "+")))
    formula.NT <- as.formula(paste("Surv(Y, event) ~", paste(covariates, collapse = "+")))
    # for fitting RCctree and RCcox, setting all left-truncation time as 0
    DATA.km <- DATA
    DATA.km$A <- 0

    cens_rate_list[mm] <- mean(1-DATA$event)
    cat("************",mm,"th simulation ***********\n")
    cat("Mean Censoring rate :",round(mean(cens_rate_list),2)*10^2,'\n')


    ##################--------- F1 and L2 for RCctree --------##################
    obj <- ctree(formula = formula.NT, data=DATA)
    pred <- predict(object = obj, newdata = TEST, type = 'prob')
    pred <- predict_ctree_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
    if(model == "tree"){
      result <- f1_macro_tree(obj, true_vars, tr_trmnd, newdata = TEST)
      RES.RR$RCctree[mm] <- result$recovery
      RES.F1$RCctree[mm] <- result$f1_macro
    }
    RES.L2$RCctree[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    print("=========RCctree done...")


    ##################--------- F1 and L2 for LTRCctree --------##################
    obj <- ltrccit(formula = formula, data = DATA)
    pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_nolog)
    if(model == "tree"){
      result <- f1_macro_tree(obj, true_vars, tr_trmnd, newdata = TEST)
      RES.RR$LTRCctree[mm] <- result$recovery
      RES.F1$LTRCctree[mm] <- result$f1_macro
    }
    RES.L2$LTRCctree[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    print("=========LTRcctree done...")


    ##################--------- F1 and L2 for LBRCctree --------##################
    obj <- ltrccit(formula = formula, data = DATA, lenbias = TRUE)
    pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_lb)
    if(model == "tree"){
      result <- f1_macro_tree(obj, true_vars, tr_trmnd, newdata = TEST)
      RES.RR$LBRCctree[mm] <- result$recovery
      RES.F1$LBRCctree[mm] <- result$f1_macro
    }
    RES.L2$LBRCctree[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    print("=========LBRcctree done...")



    ##################--------- L2 for RCcforest --------##################
    ## with default settings
    obj <- ltrccif(formula = formula, data = DATA.km, mtry=ceiling(sqrt(length(covariates))),
                   FUN = .pred_Surv_nolog, ntree = ntree, control = control, trace = FALSE)
    pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_nolog)
    RES.L2$RCcforest.D[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    print("=========RCcforest with default done...")

    # ## with mtry tuning
    # obj <- ltrccif(formula = formula, data = DATA.km, mtry=NULL,
    #                FUN = .pred_Surv_nolog, ntree = ntree, control = control, trace = FALSE)
    # pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_nolog)
    # RES.L2$RCcforest.T[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    # print("=========RCcforest with tuning done...")


    ##################--------- L2 for LTRCcforest --------##################
    ## with default settings
    obj <- ltrccif(formula = formula, data = DATA, mtry=ceiling(sqrt(length(covariates))),
                   FUN = .pred_Surv_nolog, ntree = ntree, control = control, trace = FALSE)
    pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_nolog)
    RES.L2$LTRCcforest.D[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    print("=========LTRCcforest with default done...")

    # ## with mtry tuning
    # obj <- ltrccif(formula = formula, data = DATA, mtry=NULL,
    #                FUN = .pred_Surv_nolog, ntree = ntree, control = control, trace = FALSE)
    # pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_nolog)
    # RES.L2$LTRCcforest.T[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    # print("=========LTRCcforest with tuning done...")


    ##################--------- L2 for LBRCcforest --------##################
    ## with default settings
    obj <- ltrccif(formula = formula, data = DATA, mtry=ceiling(sqrt(length(covariates))), lenbias = TRUE,
                   FUN = .pred_Surv_lb, ntree = ntree, control = control, trace = FALSE)
    pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_lb)
    RES.L2$LBRCcforest.D[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    print("=========LBRCcforest with default done...")

    # ## with mtry tuning
    # obj <- ltrccif(formula = formula, data = DATA, mtry=NULL, lenbias = TRUE,
    #                FUN = .pred_Surv_lb, ntree = ntree, control = control, trace = FALSE)
    # pred <- predictProb(object = obj, newdata = TEST, time.eval = time.uniq, FUN = .pred_Surv_lb)
    # RES.L2$LBRCcforest.T[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)
    # print("=========LBRCcforest with tuning done...")


    ##################--------- L2 for RCcox --------##################
    obj <- coxph(formula = formula.NT, data = DATA)
    pred <- survfit(obj, newdata = TEST)
    pred <- lapply(1:n, function(i) pred[i])
    pred <- predict_ctree_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
    RES.L2$RCcox[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)


    ##################--------- L2 for LTRCcox --------##################
    obj <- coxph(formula = formula, data = DATA)
    pred <- survfit(obj, newdata = TEST)
    pred <- lapply(1:n, function(i) pred[i])
    pred <- predict_ctree_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
    RES.L2$LTRCcox[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)


    ##################--------- L2 for LBRCcox --------##################
    ## EE1
    DATA.cens <- data.frame(Y=DATA$Y, event_C=1-DATA$event)
    km_C <- survfit(Surv(Y, event_C) ~ 1, data = DATA.cens)
    time_Y_m_A <- DATA$Y-DATA$A
    surv_C <- getsurv(obj = km_C, times = time_Y_m_A)
    W1 <- 1/(DATA$Y * surv_C)
    formula.lbcox <- as.formula(paste("Surv(Y, event) ~", paste(c(covariates,"offset(log(W1))"), collapse = "+")))
    obj <- coxph(formula = formula.lbcox, data = DATA, subset = event == 1)
    # different result from coxphlb(formula = formula, data=DATA)
    pred <- survfit(obj, newdata = TEST)
    pred <- lapply(1:n, function(i) pred[i])
    pred <- predict_ctree_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
    RES.L2$LBRCcox.I[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)

    ## EE2
    integrand <- function(t) getsurv(km_C,t)
    Integrate <- function(f, time.pnt){
      f.value <- sapply(time.pnt, f)
      result <- diff(time.pnt)%*%(f.value[-length(f.value)] + f.value[-1])/2
      return(result)
    }
    sort_Y <- c(0,sort(DATA$Y))
    Y_list <- sapply(1:n, function(i) sort_Y[sort_Y<=DATA$Y[i]])
    W2 <- 1/(sapply(1:n, function(i) Integrate(f = integrand, time.pnt = Y_list[[i]])))
    formula.lbcox <- as.formula(paste("Surv(Y, event) ~", paste(c(covariates,"offset(log(W2))"), collapse = "+")))
    obj <- coxph(formula = formula.lbcox, data = DATA, subset = event == 1)
    pred <- survfit(obj, newdata = TEST)
    pred <- lapply(1:n, function(i) pred[i])
    pred <- predict_ctree_cox(pred = pred, data = TEST[,c("A","Y","event","id")], time.eval = time.uniq)
    RES.L2$LBRCcox.II[mm] <- Loss.func(pred$survival.probs, TEST.true, time.uniq)


    boxplot(RES.L2$RCctree,
            RES.L2$LTRCctree,
            RES.L2$LBRCctree,
            RES.L2$RCcforest.D,
            RES.L2$LTRCcforest.D,
            RES.L2$LBRCcforest.D,
            # RES.L2$RCcforest.T,
            # RES.L2$LTRCcforest.T,
            # RES.L2$LBRCcforest.T,
            RES.L2$RCcox,
            RES.L2$LTRCcox,
            RES.L2$LBRCcox.I,
            RES.L2$LBRCcox.II,
            names = c("tree","LTtree","LBtree",
                      "forest.D","LTforest.D","LBforest.D",
                      # "RCcforest.T","LTRCcforest.T","LBRCcforest.T",
                      "cox","LTcox","LBcox.I","LBcox.II")
    )


    # printing results
    cat('\n')

    cat("  RCctree F1 :",round(sum(RES.F1$RCctree,na.rm=T)/mm,3),'\n')
    cat("LTRCctree F1 :",round(sum(RES.F1$LTRCctree,na.rm=T)/mm,3),'\n')
    cat("LBRCctree F1 :",round(sum(RES.F1$LBRCctree,na.rm=T)/mm,3),'\n')
    cat('\n')
    cat("  RCctree L2 :",round(sum(RES.L2$RCctree,na.rm=T)/mm,7)*10^5,'\n')
    cat("LTRCctree L2 :",round(sum(RES.L2$LTRCctree,na.rm=T)/mm,7)*10^5,'\n')
    cat("LBRCctree L2 :",round(sum(RES.L2$LBRCctree,na.rm=T)/mm,7)*10^5,'\n')
    cat('\n')

    cat("  RCcforest.D L2 :",round(sum(RES.L2$RCcforest.D,na.rm=T)/mm,7)*10^5,'\n')
    cat("LTRCcforest.D L2 :",round(sum(RES.L2$LTRCcforest.D,na.rm=T)/mm,7)*10^5,'\n')
    cat("LBRCcforest.D L2 :",round(sum(RES.L2$LBRCcforest.D,na.rm=T)/mm,7)*10^5,'\n')

    cat("  RCcforest.T L2 :",round(sum(RES.L2$RCcforest.T,na.rm=T)/mm,7)*10^5,'\n')
    cat("LTRCcforest.T L2 :",round(sum(RES.L2$LTRCcforest.T,na.rm=T)/mm,7)*10^5,'\n')
    cat("LBRCcforest.T L2 :",round(sum(RES.L2$LBRCcforest.T,na.rm=T)/mm,7)*10^5,'\n')
    cat('\n')

    cat("  RCcox    L2 :",round(sum(RES.L2$RCcox,na.rm=T)/mm,7)*10^5,'\n')
    cat("LTRCcox    L2 :",round(sum(RES.L2$LTRCcox,na.rm=T)/mm,7)*10^5,'\n')
    cat("LBRCcox.I  L2 :",round(sum(RES.L2$LBRCcox.I,na.rm=T)/mm,7)*10^5,'\n')
    cat("LBRCcox.II L2 :",round(sum(RES.L2$LBRCcox.II,na.rm=T)/mm,7)*10^5,'\n')

    cat('\n')
    cat('\n')
    cat('\n')
    cat('\n')


    RES <- list(RR_trees = RES.RR, F1_macro = RES.F1, L2 = RES.L2)
    RES.name <- sprintf("LBRC_DIST_%s_MODEL_%s_N%1.0f_C%1.0f",Dist,model,n,mean_cens)
    save(RES, file=RES.name)
  }
  return(RES)
}

# global setting
n <- 100
ksi <- 500
cov_set_num <- 2
M <- 100

############ tree #############
#### exponential
# exponential tree with censoring 20%
# Dist = "Exp"; model = "tree"; exp_cens_rate = 0.065; mean_cens <- 20
# exponential tree with censoring 50%
# Dist = "Exp"; model = "tree"; exp_cens_rate = 0.31; mean_cens <- 50
#### weibull increasing
# weibull increasing tree with censoring 20%
# Dist = "WI"; model = "tree"; exp_cens_rate = 0.085; mean_cens <- 20
# weibull increasing tree with censoring 50%
# Dist = "WI"; model = "tree"; exp_cens_rate = 0.33; mean_cens <- 50
#### weibull decreasing
# weibull decreasing tree with censoring 20%
# Dist = "WD"; model = "tree"; exp_cens_rate = 0.07; mean_cens <- 20
# weibull decreasing tree with censoring 50%
Dist = "WD"; model = "tree"; exp_cens_rate = 0.32; mean_cens <- 50


############ linear #############
#### exponential
# exponential linear with censoring 20%
# Dist = "Exp"; model = "linear"; exp_cens_rate = 0.002; mean_cens <- 20
# exponential linear with censoring 50%
# Dist = "Exp"; model = "linear"; exp_cens_rate = 0.018; mean_cens <- 50
#### weibull increasing
# weibull increasing linear with censoring 20%
# Dist = "WI"; model = "linear"; exp_cens_rate = 0.95; mean_cens <- 20
# weibull increasing linear with censoring 50%
# Dist = "WI"; model = "linear"; exp_cens_rate = 9; mean_cens <- 50
#### weibull decreasing
# weibull decreasing linear with censoring 20%
# Dist = "WD"; model = "linear"; exp_cens_rate = 0.24; mean_cens <- 20
# weibull decreasing linear with censoring 50%
# Dist = "WD"; model = "linear"; exp_cens_rate = 2.8; mean_cens <- 50


############ nonlinear #############
#### exponential
# exponential nonlinear with censoring 20%
# Dist = "Exp"; model = "nonlinear"; exp_cens_rate = 0.026; mean_cens <- 20
# exponential nonlinear with censoring 50%
# Dist = "Exp"; model = "nonlinear"; exp_cens_rate = 0.14; mean_cens <- 50
#### weibull increasing
# weibull increasing nonlinear with censoring 20%
# Dist = "WI"; model = "nonlinear"; exp_cens_rate = 0.2; mean_cens <- 20
# weibull increasing nonlinear with censoring 50%
# Dist = "WI"; model = "nonlinear"; exp_cens_rate = 1.1; mean_cens <- 50
#### weibull decreasing
# weibull decreasing nonlinear with censoring 20%
# Dist = "WD"; model = "nonlinear"; exp_cens_rate = 0.05; mean_cens <- 20
# weibull decreasing nonlinear with censoring 50%
# Dist = "WD"; model = "nonlinear"; exp_cens_rate = 0.33; mean_cens <- 50




tmp <- Pred_funct(n = n, Dist = Dist, model = model, cov_set_num = cov_set_num, ksi=ksi,
                  exp_cens_rate = exp_cens_rate, mean_cens = mean_cens, M = M)



