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


#' #' Use to generate simulated data where the
#' #' survival time follows a weibull distribution
#' #' 
#' #' source: https://stats.stackexchange.com/questions/135124/
#' #' how-to-create-a-toy-survival-time-to-event-data-with-right-censoring
#' #' 
#' #' @param n number of samples
#' #' @param p number of features
#' #' @param sigma strength of noise
#' #' @param beta vector of strength of covariate relationship with 
#' simulate_data_weibull <- function(n, p, sigma, beta, status){
#'   
#'   # sample features from gaussian
#'   x <- matrix(rnorm(n*p),n,p)
#'   
#'   # time linear relationship with features
#'   time <- as.vector(x %*% beta + sigma * rnorm(n))
#'   time <- time - min(time) + 1
#'   y <- cbind(time = time, status = status)
#'   
#'   # compute signal to noise ratio
#'   snr <- var(x %*% beta) / var(time - x %*% beta)
#'   
#'   return(list(x=x, y=y, snr=snr))
#' }
#' 
#' 
#' 
#' 
#' 
#' #' Use to sample without replacement from real data
#' #' where we can't generate new data
#' #' 
#' sample_data <- function(){
#'   
#' }
#' 
#' 
#' split_train_test <- function(data, ntrain){
#'   
#' }


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
                            mc.cores=mc.cores)
  
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
    outncv <- ncv_repeated(x_train,   #don't need to keep track of outncv
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
    mse_sqrt <- sqrt(mse)   # get NaN when mse_a < mse_b
    
    # store ncv output
    err_ncv[ii] <- mean(outncv$errin)
    mse_ncv[ii] <- mse
    sd_ncv[ii] <- mse_sqrt
  }
  
  # final output list
  outlist <- list(err_test=err_test, 
                  err_cv=err_cv, 
                  err_ncv=err_ncv, 
                  sd_cv=sd_cv, 
                  sd_ncv=sd_ncv, 
                  mse_ncv=mse_ncv, 
                  simulation_params=simulation_params)
  
  return(outlist)
}






source(here::here("R/ncv.R"))
registerDoMC(cores = 3)

set.seed(123)

N_SIM=10
N_TRAIN=100
N_TEST=1000
P=10
SIGMA=3
BETA = c(rep(2,4),rep(0,P-4))
STATUS_TRAIN = rep(1,N_TRAIN)
STATUS_TEST = rep(1,N_TEST)
MC_CORES = 3

data_train <- simulate_data_base_case_nsim(nsim=N_SIM, n=N_TRAIN, p=P, sigma=SIGMA, beta=BETA, status=STATUS_TRAIN)
data_test <- simulate_data_base_case_nsim(nsim=N_SIM, n=N_TEST, p=P, sigma=SIGMA, beta=BETA, status=STATUS_TEST)

simulation_test <- simulation(data_train, data_test, nsim=N_SIM, nreps=400, nfolds=10, alpha=0.10, verbose=T, mc.cores=MC_CORES)
