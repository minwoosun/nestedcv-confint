library(survival)
library(glmnet)
library(dplyr)
library(ggplot2)


generate_sim_data = function(n, p, beta, sigma){
  
  x=matrix(rnorm(n*p),n,p)
  time=as.vector(x%*%beta+sigma*rnorm(n))
  time=time-min(time)+1
  status=rep(1,n)
  y=cbind(time=time,status=status)
  
  list(x, y)
}


ntrain=100
p=10
alpha=0.1
ntests = c(50, 100, 200, 500, 1000)
beta=c(rep(2,4),rep(0,p-4))
beta_null=rep(0,p)
sigma=5
nfolds=10
mc.cores=4
nbias=500

set.seed(888)

# generate train data
data_train = generate_sim_data(n=ntrain, p=p, beta=beta, sigma=sigma)
x_train = data_train[[1]]
y_train = data_train[[2]]

# generate train data no signal
data_train_null = generate_sim_data(n=ntrain, p=p, beta=beta_null, sigma=sigma)
x_train_null = data_train_null[[1]]
y_train_null = data_train_null[[2]]

# fit model on train data
outcv = cv.glmnet(x_train,y_train,family="cox",standardize=F,keep=T)
outcv_null = cv.glmnet(x_train_null,y_train_null,family="cox",standardize=F,keep=T)

# generate and evaluate growing ntest 
err_PL = c()  # -2log(PL)
err_PL_n = c() # -2log(PL) / n
err_PL_nlogn = c() # -2log(PL) / nlog(n)
err_PL_diff = c() # -2logPL(0)-(-2logPL)
err_C = c()

for (ii in 1:length(ntests)){
  
  data_test = generate_sim_data(n=ntests[ii], p=p, beta=beta, sigma=sigma)
  xtest = data_test[[1]]
  ytest = data_test[[2]]
  
  # data_test_null = generate_sim_data(n=ntests[ii], p=p, beta=beta_null, sigma=sigma)
  # xtest_null = data_test_null[[1]]
  # ytest_null = data_test_null[[2]]
  
  
  err_PL[ii] = assess.glmnet(outcv, newx=xtest, newy=ytest,s=outcv$lambda.min)$deviance
  err_PL_n[ii] = assess.glmnet(outcv, newx=xtest, newy=ytest,s=outcv$lambda.min)$deviance/ntests[ii]
  err_PL_nlogn[ii] = assess.glmnet(outcv, newx=xtest, newy=ytest,s=outcv$lambda.min)$deviance/(ntests[ii]*log(ntests[ii]))
  err_PL_diff[ii] = (assess.glmnet(outcv_null, newx=xtest, newy=ytest,s=outcv_null$lambda.min)$deviance - assess.glmnet(outcv, newx=xtest, newy=ytest,s=outcv$lambda.min)$deviance) / ntests[ii] # (-2logPL(0)-(-2logPL)) / n
  err_C[ii] =  assess.glmnet(outcv, newx=xtest, newy=ytest,s=outcv$lambda.min)$C # c-index
  
  
}

# plotting
n_err = 3
plt.ntests = rep(ntests, n_err)
plt.err = c(err_PL_n, err_PL_nlogn, err_PL_diff, err_C)
plt.type = c(rep("err_PL_n",length(ntests)), rep("err_PL_nlogn",length(ntests)),  rep("err_PL_diff",length(ntests)),  rep("err_C",length(ntests)) )
plt.df = cbind(plt.ntests, plt.err, plt.type) %>% data.frame
plt.df$plt.ntests = plt.df$plt.ntests %>% as.numeric
plt.df$plt.err = plt.df$plt.err %>% as.numeric
names(plt.df)[3] = "error_type"


Line = ggplot(plt.df, aes(x=plt.ntests, y=plt.err, color=error_type)) +
  geom_line(size=2) +
  geom_point(size=3) + 
  xlab("Sample Size") +
  ylab("Test Error") + scale_color_manual(labels = c("log(PL)/nlog(n)", "log(PL0)-log(PL)", "log(PL)/n", "c-index"), values=c( "#ffaa00", "#0091ff", "brown1", "#6cc46c"), name = "Error type") + theme_bw() + theme(     panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.5), panel.grid.minor = element_line(colour = "grey", size = 0.5), legend.text=element_text(size=12)) # +  theme_minimal()

# panel.grid.major = element_line(colour = "#808080"),

# #6cc46c green
# #ec3231 red
# #ffaa00 yellow
# #0091ff blue
Line

# plotting
n_err = 3
plt.ntests = rep(ntests, n_err)
plt.err = c(err_PL_n, err_PL_nlogn, err_PL_diff)
plt.type = c(rep("err_PL_nlogn",length(ntests)),  rep("err_PL_diff",length(ntests)),rep("err_PL_n",length(ntests)) )
plt.df = cbind(plt.ntests, plt.err, plt.type) %>% data.frame
plt.df$plt.ntests = plt.df$plt.ntests %>% as.numeric
plt.df$plt.err = plt.df$plt.err %>% as.numeric
names(plt.df)[3] = "error_type"


Line = ggplot(plt.df, aes(x=plt.ntests, y=plt.err, color=error_type)) +
  geom_line(size=2) +
  geom_point(size=3) + 
  xlab("Sample size") +
  ylab("Log partial likelihood") + scale_color_manual(labels = c("log(PL)/nlog(n)", "log(PL0)-log(PL)", "log(PL)/n"), values=c( "#ffaa00", "#0091ff", "brown1"), name = "Error type") + theme_bw() + theme(     panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.5), panel.grid.minor = element_line(colour = "grey", size = 0.5), legend.text=element_text(size=12)) # +  theme_minimal()

# panel.grid.major = element_line(colour = "#808080"),

# #6cc46c green
# #ec3231 red
# #ffaa00 yellow
# #0091ff blue
Line
