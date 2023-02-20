
generate_sim_data <- function(n, p, beta, sigma){

  x <- matrix(rnorm(n*p),n,p)
  time <- as.vector(x%*%beta+sigma*rnorm(n))
  time <- time-min(time)+1
  status <- rep(1,n)
  y <- cbind(time=time,status=status)

  list(x, y)
}


experiment_simulation_base <- function(n, p, ntest, beta, sigma, nsim, nreps=10, nfolds=10, alpha=0.10, verbose=T, mc.cores=1){

  err = errcv = sd = sd.ncv = sd.ncva = sd.ncvb = rep(NA, nsim)
  outncv <- vector("list",nsim)

  for(ii in 1:nsim){

    if(verbose){
      print(paste("############ Sim: ", ii,"############"))
    }

    # generate train data
    data_train <- generate_sim_data(n=n, p=p, beta=beta, sigma=sigma)
    x_train <- data_train[[1]]
    y_train <- data_train[[2]]

    data_test <- generate_sim_data(n=ntest, p=p, beta=beta, sigma=sigma)
    x_test <- data_test[[1]]
    y_test <- data_test[[2]]

    # temp=coxph(formula = Surv(time, status) ~ x) #just for info on the signal strength

    fold_id <- sample(rep(1:nfolds,n/nfolds))

    outcv <- cv.glmnet(x_train,y_train,family="cox",standardize=F,foldid=fold_id,type.measure="C",keep=T)
    errcv[ii] <- max(outcv$cvm)
    err[ii] <- assess.glmnet(outcv,newx=x_test,newy=y_test,s=outcv$lambda.min)$C
    lamhat <- outcv$lambda.min
    sd[ii] <- outcv$cvsd[which.max(outcv$cvm)]

    outncv[[ii]] <- ncv(x_train,y_train,lamhat,nreps=nreps,nfolds=nfolds, mc.cores=mc.cores)
    sd.ncv[ii] <- sqrt(mean((outncv[[ii]]$errin-outncv[[ii]]$errout)^2)-mean(outncv[[ii]]$errout.var))
    sd.ncva[ii] <- mean((outncv[[ii]]$errin-outncv[[ii]]$errout)^2)
    sd.ncvb[ii] <- mean(outncv[[ii]]$errout.var)
  }

  list(err, errcv, sd, sd.ncv, outncv, sd.ncva, sd.ncvb)
}


experiment_real <- function(x, y, nsim, ntrain, nreps=10, nfolds=10, alpha=0.10, mc.cores=4, verbose=F){
  #-----------------------------------------------------
  # Function to run CV and NCV on real data

  # INPUTS
  #   x      := matrix of covariates
  #   y      := matrix of status and time
  #   ntrain := number of samples for training set
  #   nreps  := number of reps for nested CV
  #   nfolds := number of folds for CV
  #   alpha  := confidence level

  # OUTPUTS
  #   err    := vector
  #   errcv  := vector
  #   sd     := vector
  #   sd.ncv := vector
  #   outncv := list
  #-----------------------------------------------------

  err = errcv = sd = sd.ncv = sd.ncva = sd.ncvb = errnull = rep(NA, nsim)
  outncv <- vector("list",nsim)

  for(ii in 1:nsim){

    if(verbose){
      print(paste("############ Sim: ", ii,"############"))
    }

    # random train-test split
    index_train <- sample(nrow(x), ntrain)
    x_train <- x[index_train,]
    y_train <- y[index_train,]
    x_test <- x[-index_train,]
    y_test <- y[-index_train,]

    # balance folds
    y_df <- data.frame(y_train)
    y_splitted <- fold(y_df, k=nfolds, method="n_rand", cat_col="status")
    fold_id <- y_splitted$.folds %>% as.character %>% as.integer()
    # fold_id <- sample(rep(1:n_folds,nrow(x_train)/n_folds))

    # naive CV
    outcv <- cv.glmnet(x_train,y_train,family="cox",standardize=F,foldid=fold_id,type.measure="C",keep=T)
    errcv[ii] <- max(outcv$cvm)
    err[ii] <- assess.glmnet(outcv,newx=x_test,newy=y_test,s=outcv$lambda.min)$C
    lamhat <- outcv$lambda.min
    sd[ii] <- outcv$cvsd[which.max(outcv$cvm)]

    # nested CV
    outncv[[ii]] <- ncv(x_train,y_train,lamhat,nreps=nreps, nfolds=nfolds, mc.cores=mc.cores)
    sd.ncv[ii] <- sqrt(mean((outncv[[ii]]$errin-outncv[[ii]]$errout)^2)-mean(outncv[[ii]]$errout.var))
    sd.ncva[ii] <- mean((outncv[[ii]]$errin-outncv[[ii]]$errout)^2)
    sd.ncvb[ii] <- mean(outncv[[ii]]$errout.var)
  }

  list(err, errcv, sd, sd.ncv, outncv, sd.ncva, sd.ncvb)
}