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
formula.lvm <- Y~X1+X2+X3+X4
lvm.model <- lvm(list(formula.lvm))
df.data <- sim(lvm.model,n)
plvm.model <- penalize(lvm.model)

lvm.fit <- estimate(lvm.model,  data = df.data,
                      control = list(constrain = TRUE))


dataX_orth <- as.matrix(df.data[,all.vars(formula.lvm)[-1]])
ones <- rep(1,nrow(dataX_orth))
dataX_orth <- dataX_orth - cbind(ones) %*% solve(crossprod(ones), crossprod(ones, dataX_orth))
dataX_orth <- sweep(dataX_orth, MARGIN = 2, FUN = "/", 
                    STATS = sqrt(apply(dataX_orth, 2, var)*(nrow(dataX_orth)-1)/nrow(dataX_orth))
)

#### check regularization path ####
penalized.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = dataX_orth, 
                           steps = "Park")

resLassoPath <- LassoPath(dataX = dataX_orth, 
                          dataY = df.data[,all.vars(formula.lvm)[1]], 
                          iter_max = 100)

# n.knot <- nrow(resLassoPath)
# diffBeta <- unlist(lapply(penalized.fit[1:n.knot],function(x){x@lambda1})) - resLassoPath$lambda
# diffLambda <- t(data.frame(lapply(penalized.fit[1:n.knot],function(x){c(x@unpenalized, x@penalized)}))) - resLassoPath[,-1]

df.data_orth <- data.frame(Y = df.data[,all.vars(formula.lvm)[1]],
                           dataX_orth)
p1lvm.fit <- estimate(plvm.model,  data = df.data_orth, lambda1 = "Park",
                      control = list(constrain = TRUE, iter.max = 1000))


coef(penalized.fit[[1]])
penalized.fit[[3]]@lambda1
penalized.fit[[3]]@nuisance$sigma2

p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 =  572.5238,
                      control = list(constrain = TRUE, iter.max = 1000, fix.sigma = TRUE))
coef(p1lvm.fit)




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
penalized.fit <- penalized(df.data$Y,df.data[,c("X1", "X2", "X3", "X4")], 
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

