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


Pred_funct <- function(Dist="WI", n=300, exp_cens_rate=0.1, ksi=100,
                       M=500, mean_cens=20){
  # All results are stored in the results directory
  final_location <- "./results_tmp/test_unbiasedness"
  # final_location <- "./results/test_unbiasedness""

  if (!dir.exists(final_location)) {
    dir.create(final_location, recursive = TRUE, showWarnings = FALSE)
  }

  setwd(final_location)

  string = paste0("LBRC_UNBIASED_DIST_",Dist,
                  "_C",mean_cens,
                  "_sim",M)

  loaded <- tryCatch({
    load(file = string)
    F
  }, error = function(e) {
    T
  })

  if(loaded){
    M_init <- 1
    RES <- list()

    # tree Root Var
    RES$RV$LBRCctreeC <- rep(NA,M)
    RES$RV$LBRCctreeF <- rep(NA,M)

    # mean censoring rate
    RES$mean_cens_rate <- c()
  }else{
    M_init <- sum(!is.na(RES$RV$LBRCctreeF))
    if(M_init==M){
      cat(string,"finished-!! Summaries :","\n")
      cat("MCLE :", table(RES$RV$LBRCctreeC), '\n')
      cat("MFLE :", table(RES$RV$LBRCctreeF), '\n')
      return(RES)
    }
    M_init <- M_init + 1
  }

  # important!! This forces splitting
  control = partykit::ctree_control(
    mincriterion = 0,
    stump        = T,
    minsplit     = 2,
    minbucket    = 1
  )

  # tuning parameter for LBRC-CIT-(F,-)
  vardi_tune <- list(eps=1e-7,max_iter=20)

  for(mm in M_init:M){
    ## create the simulation dataset
    set.seed(mm)
    DATA <- LBRC.generate_one_sample(Dist = Dist,
                                     n = n,
                                     exp_cens_rate = exp_cens_rate,
                                     ksi = ksi)

    covariates <- setdiff(colnames(DATA), c("A","Y","event","id"))
    formula <- as.formula(paste("Surv(A, Y, event) ~", paste(covariates, collapse = "+")))

    RES$mean_cens_rate[mm] <- mean(1-DATA$event)
    cat("************",mm,"th simulation with mean censoring",round(mean(RES$mean_cens_rate),2),"***********\n")

    ##################--------- LBRC-CIT-(C,-) --------##################
    obj <- lbrccit(formula = formula, data = DATA,
                   perm_test_est = 'MCLE',
                   control = control)
    RES$RV$LBRCctreeC[mm] <- root_node_return(obj)
    if(mm>2) cat("LBRC-CIT-(C,-) :", table(RES$RV$LBRCctreeC),round(chisq.test(table(RES$RV$LBRCctreeC))$p.value,3), '\n')

    ##################--------- LBRC-CIT-(F,-) --------##################
    obj <- lbrccit(formula = formula, data = DATA,
                   perm_test_est="MFLE", perm_test_args = vardi_tune,
                   control = control)
    RES$RV$LBRCctreeF[mm] <- root_node_return(obj)
    if(mm>2) cat("LBRC-CIT-(F,-) :", table(RES$RV$LBRCctreeF),round(chisq.test(table(RES$RV$LBRCctreeF))$p.value,3), '\n')


    RES.name <- sprintf("LBRC_UNBIASED_DIST_%s_C%1.0f_sim%1.0f",
                        Dist,mean_cens,M)
    save(RES, file=RES.name)
  }
  return(RES)
}

######## simulation of unbiasedness of variable selection ########
n <- 300
ksi <- 100
M <- 10000
Dist <- "WI";  exp_cens_rates <- c(0.2,0.8);   mean_censs <- c(20,50)
# Dist <- "WD";  exp_cens_rates <- c(0.11,0.44); mean_censs <- c(20,50)
# Dist <- "Lgn"; exp_cens_rates <- c(0.09,0.33); mean_censs <- c(20,50)

for(c in 1:2){
  exp_cens_rate <- exp_cens_rates[c]
  mean_cens <- mean_censs[c]
  tmp <- Pred_funct(Dist=Dist, n=n, exp_cens_rate=exp_cens_rate,
                    M=M, mean_cens=mean_cens, ksi=ksi)
}


######## print the results ########
for(c in 1:2){
  exp_cens_rate <- exp_cens_rates[c]
  mean_cens <- mean_censs[c]
  string = paste0("LBRC_UNBIASED_DIST_",Dist,
                  "_C",mean_cens,
                  "_sim",M)
  load(file = string)
  tmpC <- table(RES$RV$LBRCctreeC)
  cat("===================================","\n")
  cat("LBRC-CIT-(C,-) :","\n")
  cat("table of selected variables : ")
  print(tmpC)
  cat("p-value :",round(chisq.test(tmpC)$p.value,4),"\n")
  cat("===================================","\n")

  tmpF <- table(RES$RV$LBRCctreeF)
  cat("===================================","\n")
  cat("LBRC-CIT-(F,-) :","\n")
  cat("table of selected variables : ")
  print(tmpF)
  cat("p-value :",round(chisq.test(tmpF)$p.value,4),"\n")
  cat("===================================","\n")
}

