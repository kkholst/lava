rm(list = ls())

#### load
library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(deSolve)
library(lava)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #### set the local path to the R files
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

source(file.path(path.lava,"tests","FCT.R"))
#### simulation
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)






##### Gold standard
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seqPark_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seqParkNorm_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))

#### Internal derivatives

source(file.path(path.lava,"R","Penalty_EPSODE.R"))
       
plvm.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2, 
                        control = list(constrain = FALSE, trace = TRUE))


## check
plvm.test <- estimate(plvm.model, data = df.data, lambda1 = plvm.EPSODE$opt$message[3,"lambda1"])


#### External derivatives
external.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = TRUE, trace = TRUE,
                            objective = objectiveO, gradient = gradientO, hessian = hessianO, 
                            control = list(constrain = FALSE))

external.EPSODE_free <- estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = FALSE, trace = TRUE,
                            objective = objectiveO, gradient = gradientO, hessian = hessianO, 
                            control = list(constrain = FALSE))

plvm.test <- estimate(plvm.model, data = df.data, lambda1 = plvm.EPSODE$opt$message[3,"lambda1"],
                      objective = objectiveO, gradient = gradientO, hessian = hessianO)
plvm.test <- estimate(plvm.model, data = df.data, lambda1 = plvm.EPSODE$opt$message[3,"lambda1"])

# external.EPSODE == plvm.EPSODE
# plvm.EPSODE_free != external.EPSODE_free


# graphical model: derivatives dSigma^-1 = Omega^-1 dOmega Omega^-1
