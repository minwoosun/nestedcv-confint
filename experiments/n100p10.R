start.time <- Sys.time()
library(dplyr)
library(glmnet)
library(survival)
library(groupdata2)
library(here)
#require(doMC)
#registerDoMC(cores = 32)
source(here::here("R/ncv.R"))
source(here::here("R/experiment.R"))

#outdir = "/scratch/users/minwoos/ncv/output/"
#jobname = "ngtp_n100_sigma10_nrep400_2"





n=100
p=10
alpha=0.1
ntest=1000
beta=c(rep(2,4),rep(0,p-4))
#beta=rep(0,p)
sigma=10
nsim=100  # repeat nested cv method as a whole nsim times to see how many times the confidence interval captures test error
nreps=400 # repeat ncv_single to get estimates for a & b
nfolds=10
mc.cores=1


result_base = experiment_simulation_base(n=n,
                                         p=p,
                                         alpha=alpha,
                                         ntest=ntest,
                                         beta=beta,
                                         sigma=sigma,
                                         nsim=nsim,
                                         nreps=nreps,
                                         nfolds=nfolds,
                                         mc.cores=mc.cores)
# assign output to respective variables
result_base[[length(result_base)+1]] = list(
                         c("1: n","2: p","3: alpha","4: ntest","5: beta","6: sigma","7: nsim","8; nreps","9: nfolds"),
                         n,
                         p,
                         alpha,
                         ntest,
                         beta,
                         sigma,
                         nsim,
                         nreps,
                         nfolds,
                         mc.cores)


# save(result_base, file=paste0(outdir,jobname,".RData"))
