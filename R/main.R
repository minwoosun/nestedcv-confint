ncv = function(x,y,lamhat,nfolds=10,nreps=5,mc.cores=4){

  raw<- parallel::mclapply(1:nreps, function(i){ncv.helper(x, y, nfolds=nfolds,lamhat=lamhat)}, mc.cores = mc.cores)
  errin=errcv0=errout=errout.var=NULL

  for(i in 1:nreps){
    errin=c(errin,raw[[i]]$errin)
    errcv0=c(errcv0,raw[[i]]$errcv0)
    errout=c(errout,raw[[i]]$errout)
    errout.var=c(errout.var,raw[[i]]$errout.var)

  }

  return(list(errin=errin,errcv0=errcv0,errout=errout,errout.var=errout.var))
}


ncv.helper = function(x,y,lamhat,nfolds=10){

  errin=errout=errout.var=rep(NA,nfolds)

  # balance folds
  y_df <- data.frame(y)
  y_splitted <- fold(y_df, k=nfolds, method="n_rand", cat_col="status")
  fold_id <- y_splitted$.folds %>% as.character %>% as.integer()

  # fold_id = sample(rep(seq(nfolds), length = nrow(x)))

  outcv0=cv.glmnet(x,y,family="cox",standardize=F,foldid=fold_id,type.measure="C",parallel=TRUE,lambda=c(lamhat,0))
  errcv0=outcv0$cvm[1]

  for(ii in 1:nfolds){
    cat(ii)
    out=which(fold_id==ii)
    xx=x[-out,]
    yy=y[-out,]
    new.foldid=fold_id[-out]
    new.foldid[new.foldid>ii]=new.foldid[new.foldid>ii]-1 #map back to 1...(nfolds-1)

    outcv=cv.glmnet(xx,yy,family="cox",standardize=F,foldid=new.foldid,type.measure="C",parallel=TRUE,lambda=c(lamhat,0))
    errin[ii]=outcv$cvm[1]

    bhat=as.vector(coef(outcv,s=outcv$lambda[1]))
    #   df=data.frame(tim=yy[,1],status=yy[,2],x=I(xx))
    xx = xx %>% data.frame
    yy = yy %>% data.frame
    df = cbind(yy, xx)
    #colnames(df)[c(1,2)] = c('status', 'time') # <- problem order will matter


    fit=coxph(Surv(time,status)~.,data=df,init=bhat,control=coxph.control(iter.max=0))
    # dfout=data.frame(x=I(x[out,]),tim=y[out,1],status=y[out,2])
    dfout = cbind(data.frame(y)[out,], data.frame(x)[out,])

    concord <- concordance(fit, newdata=dfout)
    errout[ii]=concord$concordance
    errout.var[ii]=concord$var
  }

  return(list(errin=errin,errout=errout,errout.var=errout.var,errcv0=errcv0))
}
