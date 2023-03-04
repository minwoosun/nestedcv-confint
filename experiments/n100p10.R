library(here)

source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))

# simulation parameters
SIMULATION_SEED <- 123
N_SIM <- 50
N_REP <- 100
N_TRAIN <- 100
N_TEST <- 1000
N_FOLDS <- 10
P <- 10
SIGMA <- 3
BETA <- c(rep(2, 4), rep(0, P-4))
STATUS_TRAIN <- rep(1, N_TRAIN)
STATUS_TEST <- rep(1, N_TEST)
ALPHA <- 0.10
MC_CORES <- 3
# OUTPUT_DIR

#####################################
#         Run simulation            #
#####################################

registerDoMC(cores = MC_CORES)
set.seed(SIMULATION_SEED)
#start.time <- Sys.time()

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

# save(simulation_result, file=paste0(OUTPUT_DIR,jobname,".RData"))

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

miscoverage_cv <- check_miscoverage(ci_cv, err_test)   # regular conf int shouldn't get this bad cov....
#err_test > ci_cv ["up"]
miscoverage_cnv <- check_miscoverage(ci_ncv, err_test) 

# will be NAs, have to repeat until reach desired nsim (eg. 100)
# (ran 50 iterations, only 31 not na)

# nested CV confidence intervals
