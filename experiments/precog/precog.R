here::i_am("experiments/precog/precog.R")

library(here)
library(dplyr)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))


#####################################
#           Parameters              #
#####################################

# simulation parameters
RUN_SIMULATION <- TRUE
SIMULATION_SEED <- 123
N_SIM <- 150
ALPHA <- 0.10
MC_CORES <- 32

# data & ncv parameters
N_REP <- 200
N_FOLDS <- 5
N_TRAIN <- 150

# output parameters (input_path, output_path, jobname)
args <- commandArgs(trailingOnly = TRUE)
INPUT_DIR <- args[1]
OUTPUT_DIR <- args[2]
JOBNAME <- args[3]


#######################################
#       Load and preprocess data      #
#######################################
load(here::here("data/precog_train.RData"))

###################################
#  Create NSIM train-test splits  #
###################################
split_df_list <- split_train_test_nsim_times(df_list, N_TRAIN, N_SIM)
data_train <- split_df_list[["train"]]
data_test <- split_df_list[["test"]]


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
  
  simulation_result <- simulation(data_train, 
                                  data_test, 
                                  nsim=N_SIM, 
                                  nreps=N_REP, 
                                  nfolds=N_FOLDS, 
                                  alpha=ALPHA, 
                                  verbose=TRUE, 
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
  # load(file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
  # load(here::here("data/colon.RData"))
}


#####################################
#  Output and Confidence Intervals  #
#####################################

df_result <- create_result_dataframe(simulation_result, na.rm=TRUE)
err_test <- df_result[["err_test"]]
n_sim <- nrow(df_result)

# standard CV confidence intervals
ci_cv <- confidence_interval_cv(err_cv=df_result[["err_cv"]],
                                sd_cv=df_result[["sd_cv"]],
                                alpha=ALPHA
)

ci_ncv <- confidence_interval_ncv(err_cv=df_result[["err_cv"]],
                                  err_ncv=df_result[["err_ncv"]],
                                  sd_ncv=df_result[["sd_ncv"]],
                                  alpha=ALPHA,
                                  nfolds=N_FOLDS,
                                  bias=TRUE
)

miscoverage_cv <- check_miscoverage(ci_cv, err_test)
miscoverage_ncv <- check_miscoverage(ci_ncv, err_test) 

paste0("miscoverage cv: ", miscoverage_cv)
paste0("miscoverage ncv: ", miscoverage_ncv)
