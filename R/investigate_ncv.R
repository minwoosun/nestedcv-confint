library(here)

source(here::here("R/ncv.R"))
registerDoMC(cores = 1)


generate_sim_data <- function(n, p, beta, sigma){
  
  x <- matrix(rnorm(n*p),n,p)
  time <- as.vector(x%*%beta+sigma*rnorm(n))
  time <- time-min(time)+1
  status <- rep(1,n)
  y <- cbind(time=time,status=status)
  
  list(x, y)
}


nsim=5
n=100
p=10
alpha=0.1
ntest=100
beta=c(rep(2,4),rep(0,p-4))
#beta=rep(0,p)
sigma=10
nreps=10
nfolds=10
mc.cores=1

# initialize
err = errcv = sd = sd.ncv = sd.ncva = sd.ncvb = rep(NA, nsim)
outncv <- vector("list",nsim)

for(ii in 1:nsim){
  
  cat("\n", "[ Sim: ", ii,"]", "\n")
  
  # generate train data
  data_train <- generate_sim_data(n=n, p=p, beta=beta, sigma=sigma)
  x_train <- data_train[[1]]
  y_train <- data_train[[2]]
  
  # generate test data
  data_test <- generate_sim_data(n=ntest, p=p, beta=beta, sigma=sigma)
  x_test <- data_test[[1]]
  y_test <- data_test[[2]]
  
  # temp=coxph(formula = Surv(time, status) ~ x) #just for info on the signal strength
  
  fold_id <- sample(rep(1:nfolds,n/nfolds))
  
  # run standard cross-validation
  outcv <- cv.glmnet(x_train,y_train,family="cox",standardize=F,foldid=fold_id,type.measure="C",keep=T)
  errcv[ii] <- max(outcv$cvm)
  err[ii] <- assess.glmnet(outcv,newx=x_test,newy=y_test,s=outcv$lambda.min)$C
  lamhat <- outcv$lambda.min
  sd[ii] <- outcv$cvsd[which.max(outcv$cvm)]
  
  # run nested cross-validation
  outncv[[ii]] <- ncv_repeated(x_train,y_train,lamhat,nreps=nreps,nfolds=nfolds, mc.cores=mc.cores)
  sd.ncv[ii] <- sqrt(mean((outncv[[ii]]$errin-outncv[[ii]]$errout)^2)-mean(outncv[[ii]]$errout.var))
  sd.ncva[ii] <- mean((outncv[[ii]]$errin-outncv[[ii]]$errout)^2)
  sd.ncvb[ii] <- mean(outncv[[ii]]$errout.var)
  
}


result_base <- list(err, errcv, sd, sd.ncv, outncv, sd.ncva, sd.ncvb)

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



summary.result = function(result, na.rm=T){
  
  err = result[[1]]
  errcv = result[[2]]
  sd = result[[3]]
  sd.ncv = result[[4]]
  outncv = result[[5]]
  sd.ncva = result[[6]]
  sd.ncvb = result[[7]]
  params = result[[8]]
  
  n = params[[2]]
  p = params[[3]]
  alpha = params[[4]]
  ntest = params[[5]]
  beta = params[[6]]
  sigma = params[[7]]
  nsim = nsim.total = params[[8]]
  nreps = params[[9]]
  nfolds = params[[10]]
  
  names = c('miscovL', 'miscovU', 'meanL', 'meanU', 
            'meanSE', 'meanMSEa', 'meanMSEb')
  ci.results = matrix(NA, nrow=3,ncol=7)
  
  if(na.rm==T){
    err = err[!is.na(result[[4]])]
    errcv = errcv[!is.na(result[[4]])]
    sd = sd[!is.na(result[[4]])]
    sd.ncv = sd.ncv[!is.na(result[[4]])]
    sd.ncva = sd.ncva[!is.na(result[[4]])]
    sd.ncvb =  sd.ncvb[!is.na(result[[4]])]
    nsim=length(sd)
  }
  
  
  # naive CV confidence intervals
  ci=cbind(errcv-qnorm(1-alpha/2)*sd, errcv+qnorm(1-alpha/2)*sd)
  ci.results[1,c(1,2)] = c(mean(err<ci[,1]),mean(err>ci[,2]))
  ci.results[1,c(3,4)] = colMeans(ci)
  ci.results[1,5] = mean(sd)
  ci.results[1,6] = NA
  ci.results[1,7] = NA
  
  
  # unravel nested CV point estimates
  errncv.in=matrix(NA,nsim,nreps*nfolds)
  errcv0=rep(NA,nsim)
  for(ii in 1:nsim){ 
    errncv.in[ii,]=outncv[[ii]]$errin # errin
    errcv0[ii]=mean(outncv[[ii]]$errcv0) # errout
  }
  
  # NCV point estimate
  errncv=rowMeans(errncv.in)
  
  # compute bias 
  bias = errncv-err
  biashat=((1+((nfolds-2)/nfolds)))*(errncv-errcv0)
  #biashat = errncv-errcv0
  #biashat=(errncv-errcv)
  mean(bias)
  mean(biashat)
  
  # nested CV (no bias correction)
  if (na.rm == F){
    sd.ncv.nona=sd.ncv
    sd.ncv.nona[is.na(sd.ncv)]=sd[is.na(sd.ncv)]
    sd.ncv = sd.ncv.nona
  }
  
  ci2=cbind(errncv-qnorm(1-alpha/2)*sd.ncv, errncv+qnorm(1-alpha/2)*sd.ncv)
  ci.results[2,c(1,2)] = c(mean(err<ci2[,1]),mean(err>ci2[,2]))
  ci.results[2,c(3,4)] = colMeans(ci2)
  ci.results[2,5] = mean(sd.ncv)
  ci.results[2,6] = mean(sd.ncva)
  ci.results[2,7] = mean(sd.ncvb)
  
  
  # biashat = errcv-err
  # nested CV (with bias correction)
  ci3=cbind(errncv-biashat-qnorm(1-alpha/2)*sd.ncv, errncv-biashat+qnorm(1-alpha/2)*sd.ncv)
  ci.results[3,c(1,2)] = c(mean(err<ci3[,1]),mean(err>ci3[,2]))
  ci.results[3,c(3,4)] = colMeans(ci3)
  ci.results[3,5] = mean(sd.ncv)
  ci.results[3,6] = mean(sd.ncva)
  ci.results[3,7] = mean(sd.ncvb)
  
  
  #  convert to dataframe
  ci.results.df = data.frame(ci.results)
  
  
  # first column point estimate
  point = c(mean(errcv), mean(errncv), mean(errncv-biashat))
  ci.results.df = data.frame(cbind(point,ci.results.df))
  colnames(ci.results.df) = c("meanPointEst",names)
  
  
  # summarize NCV NA sqrt(MSE)
  count_na = sum(is.na(result_base[[4]]))
  prop_na = count_na / nsim.total
  print(paste0("Sample size: ", n))
  print(paste0("Mean bias: ",mean(bias)))
  print(paste0("Mean bias-hat: ",mean(biashat)))
  print(paste0("Count of NA ncv sd: ",count_na))
  print(paste0("Proportion of NA ncv sd: ",prop_na))
  print(paste0("Mean ERR: ",mean(err)))
  
  return(ci.results.df)
}

summary.result(result_base, na.rm=T)


