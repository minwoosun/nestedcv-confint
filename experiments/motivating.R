library(here)

load(here::here("data/motivating.RData"))
#summary.result(result_base, na.rm=F)
err = result_base[[1]]
cv0 = result_base[[8]]
ncv.sd = result_base[[4]]
logx = log(seq(0.002,0.2,0.002))
alpha=0.10
ci3=cbind(rev(cv0$cvm)-qnorm(1-alpha/2)*ncv.sd, rev(cv0$cvm)+qnorm(1-alpha/2)*ncv.sd)
color.naive.ci = "#FFD100"
#color.ncv.ci = "#BEBAAE"

plot(cv0)


for(i in 1:length(err)){
  segments(logx[i],ci3[i,1],logx[i],ci3[i,2], col="grey62")
  
}
for(i in 1:length(err)){
  segments(logx[i],rev(cv0$cvup)[i],logx[i],rev(cv0$cvlo)[i], col=color.naive.ci,cex=1)
}
points(logx,ci3[,1],pch=95, col="grey62", cex=1.3)
points(logx,ci3[,2],pch=95, col="grey62", cex=1.3)
points(logx,rev(cv0$cvlo),pch=95, col=color.naive.ci, cex=1.5)
points(logx,rev(cv0$cvup),pch=95, col=color.naive.ci, cex=1.5)


points(logx,err,pch=19, col="#29D4EF", cex=0.7)
points(logx,rev(cv0$cvm),pch=19, col="red", cex=0.7)


legend(-6.2, 0.40, legend=c("True error", "CV error"),
       col=c("#29D4EF", "red"), pch=c(19,19), cex=0.8)
