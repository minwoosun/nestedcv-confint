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
#############
# load data #
#############
data(cancer, package="survival")
df = colon
df = df[complete.cases(df),]
df = df[df\$etype == 2,] # select death events
df\$rx = df\$rx  %>%  as.character
df\$rx[df\$rx == "Lev+5FU"] = 0
df\$rx[df\$rx == "Lev"] = 1
df\$rx[df\$rx == "Obs"] = 2
df\$rx = df\$rx %>% as.integer
x = df %>% select(-c(id, study, status, time, etype)) %>% data.matrix %>% as.matrix
y = df %>% select(status, time) %>% data.matrix %>% as.matrix
##################
# run CV and NCV #
##################
run=T
alpha=0.10
nsim=200
nreps=200
nfolds=10
ntrain=250
ntest = dim(x)[1] - ntrain
mc.cores=32
if (run==T){
  set.seed(888)
  exp_result = experiment_real(x=x,
                               y=y,
                               nsim=nsim,
                               ntrain=ntrain,
                               nreps=nreps,
                               nfolds=nfolds,
                               alpha=alpha,
                               mc.cores=mc.cores,
                               verbose=T)
  # assign output to respective variables
  exp_result[[length(exp_result)+1]] = list(
                           c("1: alpha ","2: ntrain", "3: ntest", "4: nsim","5: nreps","6: nfolds"),
                           alpha,
                           ntrain,
                           ntest,
                           nsim,
                           nreps,
                           nfolds)
}
save(exp_result, file=paste0(outdir,"colon3.RData"))
