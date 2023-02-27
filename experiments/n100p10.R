library(here)

source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))

# simulation parameters
SIMULATION_SEED <- 123
N_SIM=10
N_REP=10
N_TRAIN=100
N_TEST=1000
N_FOLDS=10
P=10
SIGMA=3
BETA = c(rep(2, 4), rep(0, P-4))
STATUS_TRAIN = rep(1, N_TRAIN)
STATUS_TEST = rep(1, N_TEST)
ALPHA = 0.10
MC_CORES = 3
# OUTPUT_DIR

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

####
# output and confidence intervals
####

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

miscoverage_cv <- check_miscoverage(df_ci, err_test)

# nested CV confidence intervals







#outdir = "/scratch/users/minwoos/ncv/output/"
#jobname = "ngtp_n100_sigma10_nrep400_2"



# 
# 
# n=100
# p=10
# alpha=0.1
# ntest=1000
# beta=c(rep(2,4),rep(0,p-4))
# #beta=rep(0,p)
# sigma=10
# nsim=100  # repeat nested cv method as a whole nsim times to see how many times the confidence interval captures test error
# nreps=400 # repeat ncv_single to get estimates for a & b
# nfolds=10
# mc.cores=1
# 
# 
# result_base = experiment_simulation_base(n=n,
#                                          p=p,
#                                          alpha=alpha,
#                                          ntest=ntest,
#                                          beta=beta,
#                                          sigma=sigma,
#                                          nsim=nsim,
#                                          nreps=nreps,
#                                          nfolds=nfolds,
#                                          mc.cores=mc.cores)
# # assign output to respective variables
# result_base[[length(result_base)+1]] = list(
#                          c("1: n","2: p","3: alpha","4: ntest","5: beta","6: sigma","7: nsim","8; nreps","9: nfolds"),
#                          n,
#                          p,
#                          alpha,
#                          ntest,
#                          beta,
#                          sigma,
#                          nsim,
#                          nreps,
#                          nfolds,
#                          mc.cores)
# 
# 
# # save(result_base, file=paste0(outdir,jobname,".RData"))
