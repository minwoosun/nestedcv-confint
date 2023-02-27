library(glmnet)
library(survival)
# library(groupdata2)
library(doMC)


#' Core nested cross-validation function
#' old name: ncv.helper
#' 
#' @param x data matrix
#' @param y target variable
#' @param lamhat 
#' @param nfolds number of folds for cross-validation
ncv_single = function(x, y, lamhat, nfolds=10, verbose=FALSE){
  
  # initialize output vectors
  errin = errout = errout.var = rep(NA, nfolds)

  # balance folds (balanced target)
  y_df <- data.frame(y)
  y_splitted <- fold(y_df, k=nfolds, method="n_rand", cat_col="status")
  fold_id <- as.integer(as.character(y_splitted$.folds))
  
  # run nested cross-validation
  for(ii in 1:nfolds){
    
    if(verbose){cat(paste0(ii, "."))}
    
    # define holdout fold indices
    index.out <- which(fold_id == ii)
    xx <- x[-index.out,]
    yy <- y[-index.out,]
    new.foldid <- fold_id[-index.out]
    
    # map back to 1, ..., (nfolds-1)
    new.foldid[new.foldid > ii] <- new.foldid[new.foldid > ii] - 1 
    
    # run standard k-fold cross-validation on inner folds
    outcv <- cv.glmnet(xx,
                       yy,
                       family = "cox",
                       standardize = F,
                       foldid = new.foldid,
                       type.measure = "C",
                       parallel = TRUE,
                       lambda = c(lamhat, 0)
                       )
    
    # error corresponding to lamhat model
    errin[ii] <- outcv$cvm[1]
    
    # what is the purpose of this?
    bhat <- as.vector(coef(outcv, s = outcv$lambda[1]))
    xx <- data.frame(xx)
    yy <- data.frame(yy)
    df <- cbind(yy, xx)
    #colnames(df)[c(1,2)] = c('status', 'time') # <- problem order will matter
    
    # fit a cox model using training data (replace simple lin alg representation?)
    fit <- coxph(Surv(time,status)~.,
                 data = df,
                 init = bhat,
                 control = coxph.control(iter.max = 0)
                 )
    
    # compute c-index on holdout fold using train cox model
    df_holdout_fold <- cbind(data.frame(y)[index.out,], 
                             data.frame(x)[index.out,]
                             )
    
    concord <- concordance(fit, newdata = df_holdout_fold)
    errout[ii] <- concord$concordance
    errout.var[ii] <- concord$var
  }

  outlist <- list(errin = errin,
                  errout = errout,
                  errout.var = errout.var
                  )
  
  return(outlist)
}


#' Repeated nested cross-validation function
#' Runs ncv_single() nreps times after reassigning data to new folds
#' old name: ncv
#' 
#' @param x data matrix
#' @param y target variable
#' @param lamhat 
#' @param nfolds number of folds for cross-validation
#' @param nreps number of repetitions (randomly sample x, y and run ncv)
#' @param mc.cores number of cores to use for parallel compute
ncv_repeated = function(x, y, lamhat, nfolds = 10, nreps = 5, mc.cores, verbose=FALSE){
  
  # run ncv in parallel (nreps times)
  ncv_output <- parallel::mclapply(1:nreps, 
                function(i){ncv_single(x, y, nfolds = nfolds, lamhat = lamhat, verbose)}, 
                mc.cores = mc.cores
                )
  
  # initialize output vectors
  errin = errout = errout.var = NULL
  
  # unpack output into vectors
  for(i in 1:nreps){
    errin <- c(errin, ncv_output[[i]]$errin)
    errout <- c(errout, ncv_output[[i]]$errout)
    errout.var <- c(errout.var, ncv_output[[i]]$errout.var)
  }
  
  # pack output vectors into a single list
  outlist <- list(errin = errin,
                  errout = errout,
                  errout.var = errout.var)
  
  return(outlist)
}

