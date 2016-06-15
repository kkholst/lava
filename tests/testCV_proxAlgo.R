rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(lava)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})
source(file.path(path.lava,"tests","FCT.R"))

#### 0- problematic examples ####

#### data
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

#### models
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

#### regularization path
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seqPark_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seqParkNorm_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))

#### specific knot
iterLambda <- 2
lambda_tempo <- seqParkNorm_lambda[iterLambda]

start1 <- coef(estimate(lvm.model, data = df.data))
start2 <- coef(estimate(plvm.model, lambda1 = 1e5, data = df.data))

plvm.punctual1 <- estimate(plvm.model, data = df.data, lambda1 = lambda_tempo,
                           control = list(start = start1, iter.max = 5000))

# validLVM(plvm.punctual1, penalized.PathL1[[iterLambda]])
# penalized.PathL1[[iterLambda]]@nuisance

plvm.punctual2 <- estimate(plvm.model, data = df.data, lambda1 = lambda_tempo,
                           control = list(start = start2, iter.max = 5000))

validLVM(plvm.punctual1)
validLVM(plvm.punctual2)
coef(plvm.punctual1) - coef(plvm.punctual2)



plvm.Fixed <- estimate(plvm.model, data = df.data, lambda1 = seqPark_lambda[iterLambda], fixSigma = TRUE, trace = TRUE,
                           control = list(start = start1, iter.max = 5000))
coef(penalized(Y ~  ., data = df.data, lambda1 =  seqPark_lambda[iterLambda], trace = TRUE))

