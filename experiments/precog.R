start.time <- Sys.time()
library(dplyr)
library(glmnet)
library(survival)
library(groupdata2)
require(doMC)
registerDoMC(cores = 32)
source("/home/users/minwoos/projects/ncv/code/ncv.R")
source("/home/users/minwoos/projects/ncv/code/experiment.R")
wd = "/scratch/users/minwoos/ncv/data/PRECOG_train/"
outdir = "/scratch/users/minwoos/ncv/output/"
setwd(wd)
#############
# load data #
#############
# loading annotations
temp.surv = list.files(path = wd, pattern="*.tsv")
annotations = lapply(temp.surv, read.delim)
#annotations = do.call("rbind", annotations)
#colnames(annotations) = c("sample","time", "status")
# sample time status
# loading expression matrices
temp.matrix = list.files(path = wd, pattern="*.txt|*.pcl")
for (mat in 1:length(temp.matrix)){
  filename = paste0(wd,temp.matrix[mat])
  matname = paste0("matrix",mat)
#  print(matname)
  assign(matname,read.delim(filename))
}
list.matrix = list(matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7)
print(dim(matrix1))
print("Loaded Data --------------------------")
######################
# data preprocessing #
######################
# # need to filter for rows that match between X and Y
# for (i in 1:length(list.matrix)){
#   print(length(annotations[[i]][,1]))
#   print((sum(annotations[[i]][,1] %in% colnames(list.matrix[[i]]) )) == (length(annotations[[i]][,1])) )
# 
# }
# 4th dataset has 30 samples missing in the data
#annotations[[4]] = annotations[[4]][(annotations[[4]][,1] %in% colnames(list.matrix[[4]])),] 
for (i in 1:length(annotations)){
  annotations[[i]] = annotations[[i]][(annotations[[i]][,1] %in% colnames(list.matrix[[i]])),]
}
matrices = list()
for (i in 1:length(list.matrix)){
  # remove duplicate gene names <- no way to identify among nonunique
  list.matrix[[i]] = list.matrix[[i]][!duplicated(list.matrix[[i]][,2]),]
  
  # only include samples that are in annotation
 # list.matrix[[i]] = list.matrix[[i]][, which((colnames(list.matrix[[i]]) %in% annotations[[i]][, 1]))]
  
  # transpose so rows are samples, cols are genes
  matrices[[i]] = list.matrix[[i]][,-c(1,2)] %>% t %>% data.frame
  #matrices[[i]] = list.matrix[[i]] %>% t %>% data.frame
  colnames(matrices[[i]]) = list.matrix[[i]][,2]
  
  ## need to get the order right rows of X and rows of Y
  #annotations[[i]] = annotations[[i]][order(match(annotations[[i]]\$Array, rownames(matrices[[i]]))),]
}
print(dim(matrices[[1]]))
# keep columns that intersect 
matrices2 = list()
intersecting.genes = Reduce(intersect, lapply(matrices, colnames))
for (i in 1:length(matrices)){
  # keep intersecting
  index.intersect = which(colnames(matrices[[i]]) %in% intersecting.genes)
  print(length(index.intersect))
  matrices2[[i]] = matrices[[i]][ ,index.intersect]
}
print(dim(matrices2[[1]]))
# combine all matrices into one
X = do.call("rbind", matrices2)
annotations = do.call("rbind", annotations)
#colnames(annotations) = c("sample","time", "status")
print(dim(annotations))
# remove samples not in annotation <- temp measure need to fix
probgenes = setdiff(as.character(rownames(X)),as.character(annotations\$Array))
print(length(probgenes))
print(sum(rownames(X) %in% probgenes))
print(head(probgenes))
print(head(rownames(X)))
X = X[!(rownames(X) %in% probgenes),] 
print(dim(X))
# reorder annotations to match order of X
annotations = annotations[order(match(annotations\$Array, rownames(X))),]
#setdiff(as.character(rownames(X)),as.character(annotations\$Array) )
#setdiff(as.character(annotations\$Array), as.character(rownames(X)) )
#identical(as.character(rownames(X)) ,annotations[order(match(annotations\$Array, rownames(X))),]\$Array  )
print( paste0("X and Y samples matching? ", identical(rownames(X), annotations\$Array)) )
# drop na
rows.na = which(is.na(annotations), arr.ind=TRUE)[,1]
annotations = annotations[-rows.na,]
X = X[-rows.na,]
# create response 
Y = annotations[,-1]
colnames(Y) = c('time', 'status')
#Y = Y[,c(2,1)]
# remove Y where time is <= 0
nonpos.index = which(Y[,1] <= 0) 
Y = Y[-nonpos.index,]
X = X[-nonpos.index,]
Y = data.matrix(Y)
print(paste0("X dim: ", dim(X)[1], " ", dim(X)[2]))
print(paste0("Y dim: ", dim(Y)[1], " ", dim(Y)[2]))
###########################################
# select top k features based on variance #
###########################################
X.var = apply(X, MARGIN=2, FUN=var)
names(X.var) <- NULL
X.var.df = data.frame(cbind(index=1:length(X.var), variance=X.var))
X.var.df = X.var.df %>% arrange(desc(variance))
# include FOXM1 ,KLRB1 
index.foxm1 = which(grepl("FOXM1", colnames(X), fixed = TRUE))
index.klrb1 = which(grepl("KLRB1", colnames(X), fixed = TRUE))
#top200.index = X.var.df\$index[c(1:200,index.foxm1,index.klrb1)]
top.k = 180
top.index = X.var.df\$index[1:top.k]
# normalize 
X.selected = X[,top.index]
X.selected = scale(X.selected, center=TRUE, scale=TRUE)
print("Staring CV and NCV --------------------------")
##################
# run CV and NCV #
##################
set.seed(4)
alpha=0.10
nreps=200
nfolds=5
nsim=150
mc.cores=32
ntrain=150
exp_result = experiment_real(x=X.selected, 
                             y=Y, 
                             nsim=nsim, 
                             ntrain=ntrain,
                             nreps=nreps,    #10 
                             nfolds=nfolds,  #5
                             alpha=alpha, 
                             mc.cores=mc.cores,
                             verbose=T)
# assign output to respective variables
exp_result[[length(exp_result)+1]] = list(
                       c("2: alpha ","3: ntrain", "4: nsim","5: nreps","6: nfolds"),
                       alpha,
                       ntrain,
                       nsim,
                       nreps,
                       nfolds)
save(exp_result, file=paste0(outdir,"PRECOG_result.RData"))
