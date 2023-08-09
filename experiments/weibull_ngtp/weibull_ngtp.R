here::i_am("experiments/weibull_ngtp/weibull_ngtp.R")

library(here)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))
source(here::here("R/helper.R"))

#####################################
#           Parameters              #
#####################################

# simulation parameters
RUN_SIMULATION <- TRUE
SIMULATION_SEED <- 123
N_SIM <- 150
ALPHA <- 0.10
MC_CORES <- 32

# simulation data & ncv parameters
N_REP <- 100
N_TRAIN <- 2000
N_TEST <- 5000
N_FOLDS <- 5
P <- 10
BETAS <- c(rep(0.1, 2), rep(0, P-2))
LAMBDA <- 0.01
MAXT <- 50
RATEC <- 0.02
MU <- 10
STD <- 5


# output parameters (input_path, output_path, jobname)
args <- commandArgs(trailingOnly = TRUE)
INPUT_DIR <- args[1]
OUTPUT_DIR <- args[2]
JOBNAME <- args[3]


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
  data_train <- simulate_data_exponential_nsim_times(nsim=N_SIM, 
                                                 n=N_TRAIN, 
                                                 p=P, 
                                                 beta=BETAS,
                                                 lambda=LAMBDA, 
                                                 maxt=MAXT,
                                                 rateC=RATEC,
                                                 mu=MU,
                                                 std=STD
                                                 )
  
  # simulate test data
  data_test <- simulate_data_exponential_nsim_times(nsim=N_SIM,
                                                            n=N_TEST, 
                                                            p=P, 
                                                            beta=BETAS,
                                                            lambda=LAMBDA, 
                                                            maxt=MAXT,
                                                            rateC=RATEC,
                                                            mu=MU,
                                                            std=STD
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
  save(simulation_result, file=here::here("data/weibull_ngtp.RData"))
  #save(simulation_result, file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
  #print(paste0("Saved simulation result: ",OUTPUT_DIR,JOBNAME,".RData"))

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
ci_cv1 <- confidence_interval_cv1(err_cv=df_result[["err_cv"]],
                                 sd_cv=df_result[["sd_cv"]],
                                 alpha=ALPHA
                                 )


# nested CV confidence intervals
ci_ncv_biascor <- confidence_interval_ncv(err_cv=df_result[["err_cv"]],
                                  err_ncv=df_result[["err_ncv"]],
                                  sd_ncv=df_result[["sd_ncv"]],
                                  alpha=ALPHA,
                                  nfolds=N_FOLDS,
                                  bias=TRUE
                                  )

ci_ncv_nobiascor <- confidence_interval_ncv(err_cv=df_result[["err_cv"]],
                                  err_ncv=df_result[["err_ncv"]],
                                  sd_ncv=df_result[["sd_ncv"]],
                                  alpha=ALPHA,
                                  nfolds=N_FOLDS,
                                  bias=FALSE
                                  )

miscoverage_cv1 <- check_miscoverage(ci_cv1, err_test)
miscoverage_ncv_biascor <- check_miscoverage(ci_ncv_biascor, err_test) 
miscoverage_ncv_nobiascor <- check_miscoverage(ci_ncv_nobiascor, err_test) 

paste0("miscoverage cv1: ", miscoverage_cv1)
paste0("miscoverage ncv (bias corrected): ", miscoverage_ncv_biascor)
paste0("miscoverage ncv (no bias correc): ", miscoverage_ncv_nobiascor)


plot_simulation(err_test, ci_cv1, ylim=c(0.5,0.8))
plot_simulation(err_test, ci_ncv_nobiascor , ylim=c(0.5,0.8))
plot_simulation_both(err_test, ci_cv1, ci_ncv_nobiascor , ylim=c(0.5,0.8))


index_miscov_lo <- which(err_test < ci_ncv_nobiascor[,1])


