here::i_am("experiments/base_ngtp/base_ngtp.R")

library(here)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))
source(here::here("R/helper.R"))

#####################################
#           Parameters              #
#####################################

# simulation parameters
RUN_SIMULATION <- TRUE
SIMULATION_SEED <- 888
N_SIM <- 150
ALPHA <- 0.10
MC_CORES <- 3
ACCOUNT_NAN <- TRUE
DESIRED_NSIM <- 100

# simulation data & ncv parameters
N_REP <- 100
N_TRAIN <- 1000
N_TEST <- 5000
N_FOLDS <- 5
P <- 10
BETAS <- c(rep(2, 4), rep(0, P-4))
SIGMA <- 5

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
  data_train <- simulate_data_base_case_nsim(nsim=N_SIM, 
                                                 n=N_TRAIN, 
                                                 p=P, 
                                                 beta=BETAS,
                                                 sigma=SIGMA,
                                                 status=c(rep(1, N_TRAIN))
                                                 )
  
  # simulate test data
  data_test <- simulate_data_base_case_nsim(nsim=N_SIM,
                                            n=N_TEST, 
                                            p=P, 
                                            beta=BETAS,
                                            sigma=SIGMA,
                                            status=c(rep(1, N_TEST))
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
  save(simulation_result, file=paste0(here::here("data/base_ngtp.RData")))
 # save(simulation_result, file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
#  print(paste0("Saved simulation result: ",OUTPUT_DIR,JOBNAME,".RData"))

} else {
  
  # load result
  load(file=here::here("data/base_ngtp.RData"))
}


#####################################
#  Output and Confidence Intervals  #
#####################################

if (ACCOUNT_NAN){
  # need to account for the NaN
  index.na <- which(is.na(simulation_result$sd_ncv))
  additional_rm <- (N_SIM-length(index.na)) - DESIRED_NSIM
  index_addtl_rm <- sample(c(1:N_SIM)[-index.na], size=additional_rm,  replace=FALSE)
  
  simulation_result$err_test <- simulation_result$err_test[-c(index.na, index_addtl_rm)]
  simulation_result$err_cv <- simulation_result$err_cv[-c(index.na, index_addtl_rm)]
  simulation_result$err_ncv <- simulation_result$err_ncv[-c(index.na, index_addtl_rm)]
  simulation_result$sd_cv <- simulation_result$sd_cv[-c(index.na, index_addtl_rm)]
  simulation_result$sd_ncv <- simulation_result$sd_ncv[-c(index.na, index_addtl_rm)]
  simulation_result$mse_ncv <- simulation_result$sd_mse[-c(index.na, index_addtl_rm)]
}


df_result <- create_result_dataframe(simulation_result, na.rm=TRUE)
err_test <- df_result[["err_test"]]
n_sim <- nrow(df_result)

# standard CV confidence intervals
ci_cv <- confidence_interval_cv(err_cv=df_result[["err_cv"]],
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

miscoverage_cv <- check_miscoverage(ci_cv, err_test)
miscoverage_ncv_biascor <- check_miscoverage(ci_ncv_biascor, err_test) 

paste0("miscoverage cv: ", miscoverage_cv)
paste0("miscoverage ncv (bias corrected): ", miscoverage_ncv_biascor)


plot_simulation(err_test, ci_cv, ylim=c(0.5,0.8))
plot_simulation_both(err_test, ci_cv, ci_ncv_biascor , ylim=c(0.5,0.8))


# print metrics for manuscript
mean(simulation_result$err_cv)
mean(simulation_result$err_ncv)
mean(simulation_result$sd_cv)
mean(simulation_result$sd_ncv)

