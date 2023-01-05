start.time <- Sys.time()
library(dplyr)
library(glmnet)
library(survival)
library(groupdata2)
require(doMC)
registerDoMC(cores = 32)
source("/home/users/minwoos/projects/ncv/code/ncv.R")
source("/home/users/minwoos/projects/ncv/code/experiment.R")
outdir = "/scratch/users/minwoos/ncv/output/"
jobname = "pgtn_biasfix_n100p150_nrep200_sigma10"
run=T
if (run==T){
  n=100
  p=150
  alpha=0.1
  ntest=1000
  beta=c(rep(2,4),rep(0,p-4))
  #beta=rep(0,p)
  sigma=10
  nsim=100
  nreps=200
  nfolds=10
  mc.cores=32
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
}
save(result_base, file=paste0(outdir,jobname,".RData"))

