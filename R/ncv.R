library(glmnet)
library(survival)
library(groupdata2)
library(doMC)


#' Core nested cross-validation function
#' old name: ncv.helper
#' 
#' @param x data matrix
#' @param y target variable
#' @param lamhat 
#' @param nfolds number of folds for cross-validation
ncv_single = function(x, y, lamhat, nfolds=10){
  
  # initialize output vectors
  errin = errout = errout.var = rep(NA, nfolds)

  # balance folds (balanced target)
  y_df <- data.frame(y)
  y_splitted <- fold(y_df, k=nfolds, method="n_rand", cat_col="status")
  fold_id <- as.integer(as.character(y_splitted$.folds))
  
  # run k-fold cross-validation
  outcv0 <- cv.glmnet(x, 
                      y, 
                      family = "cox", 
                      standardize = FALSE, 
                      foldid = fold_id, 
                      type.measure = "C", 
                      parallel = TRUE, 
                      lambda = c(lamhat, 0))
  
  # error corresponding to lamhat model
  errcv0 <- outcv0$cvm[1]
  
  # run nested cross- validation
  for(ii in 1:nfolds){
    
    cat(paste0(ii, "."))
    
    # define holdout fold indices
    index.out <- which(fold_id==ii)
    xx <- x[-index.out,]
    yy <- y[-index.out,]
    new.foldid <- fold_id[-index.out]
    
    # map back to 1...(nfolds-1)
    new.foldid[new.foldid>ii] <- new.foldid[new.foldid > ii] - 1 
    
    # run k-fold cross-validation
    outcv <- cv.glmnet(xx,
                       yy,
                       family="cox",
                       standardize=F,
                       foldid=new.foldid,
                       type.measure="C",
                       parallel=TRUE,
                       lambda=c(lamhat, 0))
    
    # error corresponding to lamhat model
    errin[ii] <- outcv$cvm[1]

    bhat <- as.vector(coef(outcv,s=outcv$lambda[1]))
    xx <- data.frame(xx)
    yy <- data.frame(yy)
    df <- cbind(yy, xx)
    #colnames(df)[c(1,2)] = c('status', 'time') # <- problem order will matter

    fit <- coxph(Surv(time,status)~.,data=df,init=bhat,control=coxph.control(iter.max=0))
    dfout <- cbind(data.frame(y)[index.out,], data.frame(x)[index.out,])

    concord <- concordance(fit, newdata=dfout)
    errout[ii] <- concord$concordance
    errout.var[ii] <- concord$var
  }

  outlist <- list(errin=errin,errout=errout,errout.var=errout.var,errcv0=errcv0)
  return(outlist)
}


#' Repeated nested cross-validation function
#' Runs ncv_single() nreps times
#' old name: ncv
#' 
#' @param x data matrix
#' @param y target variable
#' @param lamhat 
#' @param nfolds number of folds for cross-validation
#' @param nreps number of repetitions (randomly sample x, y and run ncv)
#' @param mc.cores number of cores to use for parallel compute
ncv_repeated = function(x, y, lamhat, nfolds = 10, nreps = 5, mc.cores){
  
  # run ncv in parallel (nreps times)
  raw <- parallel::mclapply(1:nreps, 
                            function(i){ncv_single(x, y, nfolds = nfolds, lamhat = lamhat)}, 
                            mc.cores = mc.cores)
  
  # initialize output vectors
  errin = errcv0 = errout = errout.var = NULL
  
  # unpack output into vectors
  for(i in 1:nreps){
    errin <- c(errin, raw[[i]]$errin)
    errcv0 <- c(errcv0, raw[[i]]$errcv0)
    errout <- c(errout, raw[[i]]$errout)
    errout.var <- c(errout.var, raw[[i]]$errout.var)
  }
  
  # pack output vectors into a single list
  outlist <- list(errin = errin,
                  errcv0 = errcv0,
                  errout = errout,
                  errout.var = errout.var)
  
  return(outlist)
}

