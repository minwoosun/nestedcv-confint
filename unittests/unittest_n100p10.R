here::i_am("unittests/unittest_n100p10.R")

library(here)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))


#####################################
#           Parameters              #
#####################################

# simulation parameters
RUN_SIMULATION <- TRUE
SIMULATION_SEED <- 123
N_SIM <- 10
ALPHA <- 0.10
MC_CORES <- 1

# simulation data & ncv parameters
N_REP <- 20
N_TRAIN <- 100
N_TEST <- 500
N_FOLDS <- 10
P <- 10
SIGMA <- 20
BETA <- c(rep(2, 4), rep(0, P-4))
STATUS_TRAIN <- rep(1, N_TRAIN)
STATUS_TEST <- rep(1, N_TEST)

# output parameters (input_path, output_path, jobname)
OUTPUT_DIR <- here::here("data/")
JOBNAME <- "unitTest"


#####################################
#         Run simulation            #
#####################################
# Set RUN_SIMULATION to FALSE if already have simulation results,
# and instead load the simulation results.
# When a < b for MSE estimation, you will get NaN,
# make NSIM larger to get desired NSIM as NaN will be discarded.

if (RUN_SIMULATION) {
  
  # registerDoMC(cores = MC_CORES)
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
  simulation_result[["simulation_params"]][["compute.time"]] <- compute.time
  
  # save result
 #  save(simulation_result, file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
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
ci_cv <- confidence_interval_cv1(err_cv=df_result[["err_cv"]],
                                 sd_cv=df_result[["sd_cv"]],
                                 alpha=ALPHA
                                 )

# nested CV confidence intervals
ci_ncv <- confidence_interval_ncv(err_cv=df_result[["err_cv"]],
                                  err_ncv=df_result[["err_ncv"]],
                                  sd_ncv=df_result[["sd_ncv"]],
                                  alpha=ALPHA,
                                  nfolds=N_FOLDS,
                                  bias=TRUE
                                  )

miscoverage_cv <- check_miscoverage(ci_cv, err_test)
miscoverage_ncv <- check_miscoverage(ci_ncv, err_test) 
# err_test > ci_cv ["up"]

paste0("miscoverage cv: ", miscoverage_cv)
paste0("miscoverage ncv: ", miscoverage_ncv)

unit_test_result <- (miscoverage_cv[1] == 0.125) & (miscoverage_ncv[1] == 0.125)
if(unit_test_result){
  print("PASS")
} else {
  print("FAIL")
}





