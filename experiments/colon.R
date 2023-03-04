here::i_am("experiments/colon.R")

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

# data & ncv parameters
N_REP <- 20
N_FOLDS <- 10
N_TRAIN <- 250

# output parameters
JOBNAME <- "colon"
OUTPUT_DIR <- here::here("data/simulation_output/") # <- change this for sherlock


#######################################
#       Load and preprocess data      #
#######################################
data(cancer, package="survival")
df <- colon

# data preprocessing
df <- df[complete.cases(df),]
df <- df[df$etype == 2,] # select death events
df$rx <- as.character(df$rx) 
df$rx[df$rx == "Lev+5FU"] <- 0
df$rx[df$rx == "Lev"] <- 1
df$rx[df$rx == "Obs"] <- 2
df$rx <- as.integer(df$rx)

x <- df %>% 
      select(-c(id, study, status, time, etype)) %>% 
      data.matrix %>% 
      as.matrix

y <- df %>% 
      select(status, time) %>% 
      data.matrix %>% 
      as.matrix


###################################
#  Create NSIM train-test splits  #
###################################
df_list <- list("x"=x, "y"=y)
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
  load(file=paste0(OUTPUT_DIR,JOBNAME,".RData"))
}



