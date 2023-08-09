here::i_am("experiments/colon/colon_results.R")

library(here)
source(here::here("R/ncv.R"))
source(here::here("R/simulation.R"))
source(here::here("R/helper.R"))


#####################################
#          Parameters               #
#####################################

# load data
load(here::here("data/colon.RData"))

# simulation parameters
SIMULATION_SEED <- 123
ALPHA <- 0.10
SELECT_SIZE <- TRUE
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

if (SELECT_SIZE){
  # need to account for the NaN
  index.include <- 101:(DESIRED_NSIM+100) # get second half 100

  simulation_result$err_test <- simulation_result$err_test[index.include]
  simulation_result$err_cv <- simulation_result$err_cv[index.include]
  simulation_result$err_ncv <- simulation_result$err_ncv[index.include]
  simulation_result$sd_cv <- simulation_result$sd_cv[index.include]
  simulation_result$sd_ncv <- simulation_result$sd_ncv[index.include]
  simulation_result$mse_ncv <- simulation_result$sd_mse[index.include]
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

