here::i_am("experiments/base_ngtp/base_ngtp_results.R")

library(here)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))
source(here::here("R/helper.R"))


#####################################
#           Parameters              #
#####################################

# load result
load(here::here("data/base_ngtp.RData"))

# simulation parameters
SIMULATION_SEED <- 888
ALPHA <- 0.10
ACCOUNT_NAN <- TRUE
DESIRED_NSIM <- 100
N_SIM <- simulation_result$simulation_params$nsim
N_REP <- simulation_result$simulation_params$nreps
N_TRAIN <- simulation_result$simulation_params$n_train
N_TEST <- simulation_result$simulation_params$n_test
N_FOLDS <-simulation_result$simulation_params$nfolds
P <- simulation_result$simulation_params$p

set.seed(SIMULATION_SEED)

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
plot_simulation(err_test, ci_ncv_biascor , ylim=c(0.5,0.8))
plot_simulation_both(err_test, ci_cv, ci_ncv_biascor , ylim=c(0.5,0.8))


# print metrics for manuscript
mean(simulation_result$err_cv)
mean(simulation_result$err_ncv)
mean(simulation_result$sd_cv)
mean(simulation_result$sd_ncv)

