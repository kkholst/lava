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
lambda.lvm <- 0.5
p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = lambda.lvm,
                      control = list(constrain = TRUE, iter.max = 1000))
p1lvm.fit
coef(p1lvm.fit)
lambda.penalized <- lambda.lvm*n*coef(p1lvm.fit)["y,y"]
penalized.fit <- penalized(df.data$y,df.data[,c("X1", "X2", "X3", "X4")], 
                 lambda1 = lambda.penalized)
coef(p1lvm.fit) - c(penalized.fit@unpenalized, penalized.fit@penalized,penalized.fit@nuisance$sigma2)
lambda.lvm-lambda.penalized/(n*penalized.fit@nuisance$sigma2)

require('penalized')
#### check the whole path
path.fit <- penalized(df.data$y, df.data[,c("X1", "X2", "X3", "X4")], 
                      steps = "Park" )
seq_lambda.penalized <- unlist(lapply(path.fit, function(x){x@lambda1}))
seq_sigma.penalized <- unlist(lapply(path.fit, function(x){x@nuisance$sigma2}))
seq_lambda.lvm <- seq_lambda.penalized / (n * seq_sigma.penalized)

iter_path <- 5
p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = seq_lambda.lvm[iter_path], regularizationPath = TRUE,
                      control = list(constrain = TRUE, iter.max = 1000))
coef(p1lvm.fit)
c(path.fit[[iter_path]]@unpenalized, path.fit[[iter_path]]@penalized, path.fit[[iter_path]]@nuisance$sigma2)

#### compare penalized with lava when fixing sigma
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

