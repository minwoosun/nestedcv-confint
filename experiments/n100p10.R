here::i_am("experiments/n100p10.R")

library(here)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))


#####################################
#           Parameters              #
#####################################

# simulation parameters
RUN_SIMULATION <- TRUE
SIMULATION_SEED <- 123
N_SIM <- 3
ALPHA <- 0.10
MC_CORES <- 3

# simulation data & ncv parameters
N_REP <- 100
N_TRAIN <- 100
N_TEST <- 1000
N_FOLDS <- 10
P <- 10
SIGMA <- 3
BETA <- c(rep(2, 4), rep(0, P-4))
STATUS_TRAIN <- rep(1, N_TRAIN)
STATUS_TEST <- rep(1, N_TEST)

# output parameters
JOBNAME <- "base_n100p10"
OUTPUT_DIR <- here::here("data/simulation_output/") # <- change this for sherlock


#####################################
#         Run simulation            #
#####################################
# Set RUN_SIMULATION to FALSE if already have simulation results,
# and instead load the simulation results.
# When a < b for MSE estimation, you will get NaN,
# make NSIM larger to get desired NSIM as NaN will be discarded.

if (RUN_SIMULATION) {
  
  registerDoMC(cores = MC_CORES)
  set.seed(SIMULATION_SEED)
  start.time <- Sys.time()
  
  # simulate train data (build confidence intervals)
  data_train <- simulate_data_base_case_nsim(nsim=N_SIM,
                                             n=N_TRAIN, 
                                             p=P, 
                                             sigma=SIGMA, 
                                             beta=BETA, 
                                             status=STATUS_TRAIN
                                             )
  
  # check signal-to-noise ratio
  snr_train <- check_snr(data_train, N_SIM)
  snr_train_avg <- mean(snr_train)
  
  # simulate test data
  data_test <- simulate_data_base_case_nsim(nsim=N_SIM, 
                                            n=N_TEST, 
                                            p=P, 
                                            sigma=SIGMA, 
                                            beta=BETA, 
                                            status=STATUS_TEST
                                            )
  
  # run simulations
  simulation_result <- simulation(data_train, 
                                  data_test, 
                                  nsim=N_SIM, 
                                  nreps=N_REP, 
                                  nfolds=N_FOLDS, 
                                  alpha=ALPHA, 
                                  verbose=T, 
                                  mc.cores=MC_CORES
                                  )
  
  # include simulation compute time in output
  end.time <- Sys.time()
  compute.time <- end.time - start.time
  simulation_result$simulation_params$compute.time <- compute.time
  
  # save result
  save(simulation_result, file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
  print(paste0("Saved simulation result: ",OUTPUT_DIR,JOBNAME,".RData"))

} else {
  
  # load result
  load(file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
  
}


#####################################
#  Output and Confidence Intervals  #
#####################################

df_result <- create_result_dataframe(simulation_result, na.rm=TRUE)
err_test <- df_result[["err_test"]]
n_sim <- nrow(df_result)

# standard CV confidence intervals
ci_cv <- confidence_interval_cv(err_cv=df_result["err_cv"],
                                sd_cv=df_result["sd_cv"],
                                alpha=ALPHA,
                                nsim=n_sim,
                                N_FOLDS
                                )

ci_ncv <- confidence_interval_cv(err_cv=df_result["err_ncv"],
                                sd_cv=df_result["sd_ncv"],
                                alpha=ALPHA,
                                nsim=n_sim,
                                N_FOLDS
                                )

miscoverage_cv <- check_miscoverage(ci_cv, err_test)
miscoverage_ncv <- check_miscoverage(ci_ncv, err_test) 
# err_test > ci_cv ["up"]

# will be NAs, have to repeat until reach desired nsim (eg. 100)
# (ran 50 iterations, only 31 not na)

# nested CV confidence intervals
