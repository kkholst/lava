rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})


set.seed(10)
n <- 500
lvm.model <- lvm(list(y~X1+X2+X3+X4))
df.data <- sim(lvm.model,n)
plvm.model <- penalize(lvm.model)

lvm.fit <- estimate(lvm.model,  data = df.data,
                      control = list(constrain = TRUE))

#### look at sigma
penalized.fit <- penalized(df.data$y,df.data[,c("X1", "X2", "X3", "X4")], 
                           steps = "Park")
seq_lambda <- unlist(lapply(penalized.fit,function(x){x@lambda1}))

for(iter_lambda in 1:length(seq_lambda)){

  lambda.lvm <- seq_lambda[iter_lambda]/n
  p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = lambda.lvm, 
                        control = list(constrain = TRUE, iter.max = 1000, fix.sigma = TRUE))
  
  cat("lambda = ",lambda.lvm,"\n")
  coefPenalized <- c(penalized.fit[[iter_lambda]]@unpenalized,penalized.fit[[iter_lambda]]@penalized)[c("(Intercept)","X1","X2","X3","X4")]
  print(cbind(coefLVM = coef(p1lvm.fit)[1:5],
              coefP = coefPenalized,
              diff = coef(p1lvm.fit)[1:5] - coefPenalized)
  )
  
}

