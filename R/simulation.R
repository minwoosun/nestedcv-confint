here::i_am("R/simulation.R")

library(here)
source(here::here("R/ncv.R"))


####################################################################
#                    Data generating functions                     #
####################################################################


#' Use to generate simple simulated data.
#' 
#' @param n number of samples
#' @param p number of features
#' @param sigma strength of noise
#' @param beta vector of strength of covariate relationship with time
#' @param status vector of event (e.g. 1=death, 0=censored)
simulate_data_base_case <- function(n, p, sigma, beta, status){
  
  # sample features from gaussian
  x <- matrix(rnorm(n*p),n,p)
  
  # time linear relationship with features
  time <- as.vector(x %*% beta + sigma * rnorm(n))
  time <- time - min(time) + 1
  y <- cbind(time = time, status = status)
  
  # compute signal to noise ratio
  snr <- var(x %*% beta) / var(time - x %*% beta)
  
  return(list(x=x, y=y, snr=snr))
}


#' Use to generate simple simulated data nsim times
#' 
#' @param nsim number of simulations
#' @param n number of samples
#' @param p number of features
#' @param sigma strength of noise
#' @param beta vector of strength of covariate relationship with time
#' @param status vector of event (e.g. 1=death, 0=censored)
simulate_data_base_case_nsim <- function(nsim, n, p, sigma, beta, status){
  
  simulated_data <- list()
  
  for (i in 1:nsim){
    simulated_data[[i]] <- simulate_data_base_case(n, p, sigma, beta, status)
  }
  
  return(simulated_data)
}


#' Get signal-to-noise ratio from simulated data.
#' Should expect n_sim SNR values
#' 
#' @param data list of data (x, y, snr) from 
#'             simulate_data_base_case_nsim()
#' @param nsim number of simulations
#' 
check_snr <- function(data, n_sim){
  
  snr <- c()
  for (i in 1:n_sim){
    snr[i] <- unlist(data[[i]]["snr"])
  }
  
  return(snr)
}


#' Use to generate simulated data where the survival time follows a exponential
#' distribution and the features are sample from gaussian.
#'
#' @param n number of samples
#' @param p number of features
#' @param beta vector of strength of covariate relationship with time
#' @param lambda scale parameter for exponential distribution
#' @param rateC rate parameter for Exponential distribution (use for censor)
#' e.g. simulate_data_exponential(n=100, p=10, beta=c(2,2,rep(0,8)), 
#'                                         lambda=0.01, rateC=0.001)   
simulate_data_exponential <- function(n, p, beta, lambda, rateC){
  
  # sample features from gaussian
  x <- matrix(rnorm(n*p),n,p)
  
  # sample survival time from weibull
  u <- runif(n)
  time <- -log(u) / (lambda * exp(x %*% beta))
  
  # incorporate censoring
  C <- rexp(n, rateC)
  time <- pmin(time, C)
  status <- as.numeric(time < C)
  
  y <- cbind(time = time, status = status)
  colnames(y) <- c("time", "status")
  
  return(list(x=x, y=y))
}


#' Use to generate simulated data where the survival time follows a weibull 
#' distribution and the features are binary.
#'
#' @param n number of samples
#' @param p number of features
#' @param beta vector of strength of covariate relationship with time
#' @param lambda scale parameter for Weibull distribution
#' @param rho shape parameter for Weibull distribution
#' @param rateC rate parameter for Exponential distriubtion (use for censor)
#' e.g. simulate_data_weibull(n=100, p=10, beta=c(2,2,rep(0,8)), 
#'                                         lambda=0.01, rho=0.5, rateC=0.001)   
simulate_data_weibull <- function(n, p, beta, lambda, rho, rateC){
  
  # sample p features by n bernoulli trials p times
  x <- matrix(NA, nrow=n, ncol=p)
  
  for (j in 1:p){
    x[,j] <- sample(x=c(0, 1), size=n, replace=TRUE, prob=c(0.5, 0.5))
  }
  
  # sample survival time from weibull
  u <- runif(n)
  time <- (-log(u) / (lambda * exp(x %*% beta)))^(1/rho)
  
  # incorporate censoring
  C <- rexp(n, rateC)
  time <- pmin(time, C)
  status <- as.numeric(time < C)
  
  y <- cbind(time = time, status = status)
  colnames(y) <- c("time", "status")
  
  return(list(x=x, y=y))
}


#' Use to generate simulated data where the survival time follows a Gompertz
#' distribution and the features are from gaussian.
#'
#' @param n number of samples
#' @param p number of features
#' @param beta vector of strength of covariate relationship with time
#' @param lambda scale parameter for Gompertz distribution
#' @param alpha shape parameter for Gompertz distribution
#' @param rateC rate parameter for Exponential distriubtion (use for censor)
#' e.g. simulate_data_gompertz(n=100, p=10, beta=c(2,2,rep(0,8)), 
#'                                         lambda=0.01, alpha=0.0001, rateC=0.001)   
simulate_data_gompertz <- function(n, p, beta, lambda, alpha, rateC){
  
  # sample features from gaussian
  x <- matrix(rnorm(n*p),n,p)
  
  # sample survival time from weibull
  u <- runif(n)
  time <- (1/alpha) * log(1 - (alpha * log(u)) / (lambda * exp(x %*% beta)))
  
  # incorporate censoring
  C <- rexp(n, rateC)
  time <- pmin(time, C)
  status <- as.numeric(time < C)
  
  y <- cbind(time = time, status = status)
  colnames(y) <- c("time", "status")
  
  return(list(x=x, y=y))
}


#' Use to generate simulated data where the survival time follows a exponential
#' distribution and the features are from gaussian NSIM times
#'
#' @param nsim number of simulations (i.e. repeat generating data)
#' @param n number of samples
#' @param p number of features
#' @param beta vector of strength of covariate relationship with time
#' @param lambda scale parameter for exponential distribution
#' @param rateC rate parameter for Exponential distriubtion (use for censor)
#' e.g. simulate_data_weibull_nsim_times(3, n=100, p=10, beta=c(2,2,rep(0,8)),  
#'                                        lambda=0.01, rho=0.5, rateC=0.001)
simulate_data_exponential_nsim_times <- function(nsim, n, p, beta, lambda, rateC){
  
  simulated_data <- list()
  
  for (i in 1:nsim){
    simulated_data[[i]] <- simulate_data_exponential(n, p, beta, lambda, rateC)
  }
  
  return(simulated_data) 
}


#' Use to generate simulated data where the survival time follows a weibull 
#' distribution and the features are binary NSIM times
#'
#' @param nsim number of simulations (i.e. repeat generating data)
#' @param n number of samples
#' @param p number of features
#' @param beta vector of strength of covariate relationship with time
#' @param lambda scale parameter for Weibull distribution
#' @param rho shape parameter for Weibull distribution
#' @param rateC rate parameter for Exponential distriubtion (use for censor)
#' e.g. simulate_data_weibull_nsim_times(3, n=100, p=10, beta=c(2,2,rep(0,8)),  
#'                                            lambda=0.01, rho=0.5, rateC=0.001)
simulate_data_weibull_nsim_times <- function(nsim, n, p, beta, lambda, rho, rateC){
  
  simulated_data <- list()
  
  for (i in 1:nsim){
    simulated_data[[i]] <- simulate_data_weibull(n, p, beta, lambda, rho, rateC)
  }
  
  return(simulated_data) 
}


#' Use to generate simulated data where the survival time follows a Gompertz
#' distribution and the features are from gaussian nsim times
#'
#' @param nsim number of simulations (i.e. repeat generating data)
#' @param n number of samples
#' @param p number of features
#' @param beta vector of strength of covariate relationship with time
#' @param lambda scale parameter for Gompertz distribution
#' @param alpha shape parameter for Gompertz distribution
#' @param rateC rate parameter for Exponential distriubtion (use for censor)
#' e.g. simulate_data_gompertz(n=100, p=10, beta=c(2,2,rep(0,8)), 
#'                                         lambda=0.01, alpha=0.0001, rateC=0.001)   
simulate_data_gompertz <- function(nsim, n, p, beta, lambda, alpha, rateC){
  
  simulated_data <- list()
  
  for (i in 1:nsim){
    simulated_data[[i]] <- simulate_data_gompertz(n, p, beta, lambda, alpha, rateC)
  }
  
  return(simulated_data) 
}


#' Use to sample without replacement from real data
#' where we can't generate new data
#'
#' @param data list of x and y
#' @param ntrain number of samples in the train set
split_train_test <- function(data, ntrain){
  
  x <- data[["x"]]
  y <- data[["y"]]
  
  # random sample
  index_train = sample(nrow(x), ntrain)
  
  # train-test split
  x_train = x[index_train,]
  y_train = y[index_train,]
  x_test = x[-index_train,]
  y_test = y[-index_train,]
  
  train <- list("x_train"=x_train, "y_train"=y_train)
  test <- list("x_test"=x_test, "y_test"=y_test)
  output <- list("train"=train, "test"=test)
  
  return(output)
}


#' Run split_train_test() nsim times for real data
#' where we can't generate new data
#'
#' @param data list of x and y
#' @param ntrain number of samples in the train set
#' @param nsim number of simulations
split_train_test_nsim_times <- function(data, ntrain, nsim){
  
  data_train <- list()
  data_test <- list()
  
  for (i in 1:nsim){
    
    split_data <- split_train_test(data, ntrain)
    x_train <- split_data[["train"]][["x_train"]]
    y_train <- split_data[["train"]][["y_train"]]
    x_test <- split_data[["test"]][["x_test"]]
    y_test <- split_data[["test"]][["y_test"]]
    
    data_train[[i]] <- list("x"=x_train, "y"=y_train)
    data_test[[i]] <- list("x"=x_test, "y"=y_test)
  }
  
  return(list("train"=data_train, "test"=data_test))
}


##########################################################
#                Simulation functions                    #
##########################################################


#' Simulates running vanilla cross-validation and
#' nested cross validation on simulated input data.
#' Runs nsim times to give nsim MSE estimates
#' 
#' err is c-index
#' @param data_train list of train data (x, y, snr)
#' @param data_test list of test data (x, y, snr)
#' @param nsim number of simulations
#' @param nreps number of repetition for ncv
#' @param nfolds number of folds
#' @param alpha desired miscoverage (1-alpha is desired conf int)
#' @param mc.cores number of cpu for parallel computation of ncv
#' @param verbose print simulation iteration
#' 
simulation <- function(data_train, data_test, nsim, nreps=10, nfolds=10, alpha=0.10, verbose=T, mc.cores=1){
  
  # simulation parameters
  n_train <- nrow(data_train[[1]][["x"]])
  n_test <- nrow(data_test[[1]][["x"]])
  p <- ncol(data_train[[1]][["x"]])
  simulation_params <- list(n_train=n_train,
                            n_test=n_test,
                            p=p,
                            nsim=nsim,
                            nreps=nreps,
                            nfolds=nfolds,
                            mc.cores=mc.cores,
                            compute.time=NA)
  
  # initialize outputs
  err_test = err_cv = err_ncv = sd_cv = sd_ncv = mse_ncv = rep(NA, nsim)
  
  for(ii in 1:nsim){
    
    if(verbose){print(paste("[Sim: ", ii,"]"))}
    
    # assign train and test data
    x_train <- data_train[[ii]][["x"]]
    y_train <- data_train[[ii]][["y"]]
    x_test <- data_test[[ii]][["x"]]
    y_test <- data_test[[ii]][["y"]]
    
    # generate folds 
    fold_id <- sample(rep(1:nfolds, n_train / nfolds))
    
    #####################################
    #   run standard cross-validation   #
    #####################################
    outcv <- cv.glmnet(x_train, 
                       y_train,
                       family="cox",
                       standardize=F,
                       foldid=fold_id,
                       type.measure="C",
                       keep=T
    )
    
    err_cv[ii] <- max(outcv$cvm)
    sd_cv[ii] <- outcv$cvsd[which.max(outcv$cvm)]
    lamhat <- outcv$lambda.min
    
    # compute test error
    err_test[ii] <- assess.glmnet(outcv,
                                  newx=x_test,
                                  newy=y_test,
                                  s=outcv$lambda.min)$C
    
    ###################################
    #   run nested cross-validation   #
    ###################################
    outncv <- ncv_repeated(x_train,
                           y_train,
                           lamhat, # from standard cv
                           nreps=nreps,
                           nfolds=nfolds, 
                           mc.cores=mc.cores
    )
    
    # compute MSE estimate
    mse_a <- mean((outncv$errin - outncv$errout)^2)
    mse_b <- mean(outncv$errout.var)
    mse <- mse_a - mse_b
    mse_sqrt <- sqrt(mse) # get NaN when mse_a < mse_b
    
    # store ncv output
    err_ncv[ii] <- mean(outncv$errin)
    mse_ncv[ii] <- mse
    sd_ncv[ii] <- mse_sqrt
  }
  
  # final output list
  outlist <- list(err_test = as.vector(unlist(err_test)), 
                  err_cv = as.vector(unlist(err_cv)), 
                  err_ncv = as.vector(unlist(err_ncv)), 
                  sd_cv = as.vector(unlist(sd_cv)), 
                  sd_ncv = as.vector(unlist(sd_ncv)), 
                  mse_ncv = as.vector(unlist(mse_ncv)), 
                  simulation_params = simulation_params)
  
  return(outlist)
}


#' Take output from simulate() and organize into a dataframe
#' 
#' @param result list of results from simulate()
#' @param na.rm TRUE remove rows with NaN
create_result_dataframe <- function(result, na.rm){
  
  err_test <- result[["err_test"]]
  err_cv <- result[["err_cv"]]
  err_ncv <- result[["err_ncv"]]
  sd_cv <- result[["sd_cv"]]
  sd_ncv <- result[["sd_ncv"]]
  
  df <- data.frame(cbind(err_test, err_cv, err_ncv, sd_cv, sd_ncv))
  colnames(df) <- c("err_test", "err_cv", "err_ncv", "sd_cv", "sd_ncv")
  
  if (na.rm){
    df <- df[complete.cases(df),]
  }
  
  return(df)
}


#' Compute confidence intervals for standard cross-validation output
#' 
#' @param err_cv
#' @param sd_cv
#' @param alpha desired miscoverage (1-alpha is desired conf int)
#' @param nsim
#' @param nfold number of folds
#' @param na.rm TRUE remove rows with NaN
confidence_interval_cv <- function(err_cv, sd_cv, alpha, nsim, nfold){
  z <- qnorm(1-alpha/2)
  df <- data.frame(cbind(err_cv-z*sd_cv, err_cv+z*sd_cv))
  colnames(df) <- c("lo", "up")
  return(df)
}


#' Compute confidence intervals for standard cross-validation output
#' 
#' @param confidence_intervals
#' @param test_err
check_miscoverage <- function(confidence_intervals, test_err){
  
  lower_bound_miscoverage <- mean(test_err < confidence_intervals["lo"])
  upper_bound_miscoverage <- mean(test_err > confidence_intervals["up"])
  miscoverage <- c(lower_bound_miscoverage, upper_bound_miscoverage)
  
  return(miscoverage)
}
