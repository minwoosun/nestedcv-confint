here::i_am("R/helper.R")


#' @param testError a vector of test errors
#' @param confInverval a matrix of confidence interval lower and upper bound
plot_simulation <- function(testError, confInterval, ylim=c(0.4, 1)){
  
  # plot test error and confidence intervals
  plot(testError, ylim = ylim)
  for (i in 1:length(testError)) {
    segments(i, confInterval[i, 1], i, confInterval[i, 2])
  }
  
  # color upper and lower bound
  points(confInterval[, 1], pch="_")
  points(confInterval[, 2], pch="_")
  
  # mark miscoverage with X
  index_miscov_lo <- which(testError < confInterval[,1])
  index_miscov_hi <- which(testError > confInterval[,2])
  
  points(index_miscov_lo, testError[index_miscov_lo], col="cyan3", pch=4, cex=3)
  points(index_miscov_hi, testError[index_miscov_hi], col="red", pch=4, cex=3)
  
}


plot_simulation_both <- function(testError, confInterval_cv, confInterval_ncv, ylim=c(0.4, 1)){
  
  # plot test error and confidence intervals
  plot(testError, ylim = ylim)
  for (i in 1:length(testError)) {
    segments(i, confInterval_ncv[i, 1], i, confInterval_ncv[i, 2], col="blue")
    segments(i, confInterval_cv[i, 1], i, confInterval_cv[i, 2], col="red")
  }
  
  # color upper and lower bound
  points(confInterval_ncv[, 1], pch="_", col="blue")
  points(confInterval_ncv[, 2], pch="_", col="blue")
  
  points(confInterval_cv[, 1], pch="_", col="red")
  points(confInterval_cv[, 2], pch="_", col="red")
  
}


plot_train_test_hist <- function(data_train, data_test, idx=1){
  
  # plot_train_test_hist(data_train, data_test)
  
  censor_proportion_train <- sum(data_train[[idx]]$y[,2]) / length(data_train[[idx]]$y[,2])
  censor_proportion_test <- sum(data_test[[idx]]$y[,2]) / length(data_test[[idx]]$y[,2])
  
  print("Proportion of censoring in train and test data: ")
  print(censor_proportion_train)
  print(censor_proportion_test)
  
  par(mfrow=c(2,1))
  
  hist(data_train[[idx]]$y[,1], 
       main="Historgram of Surival Time - Train",
       xlab="Survival Time")
  
  hist(data_test[[idx]]$y[,1], 
       main="Historgram of Surival Time - Test",
       xlab="Survival Time")

}
  