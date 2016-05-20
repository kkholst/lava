rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(deSolve)
library(lava)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

objectiveO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  lv <-  - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(epsilon) %*% epsilon
  
  return(-as.numeric(lv))
}

gradientO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * t(epsilon) %*% epsilon    
  gradient_beta <- + 1/(sigma2) * t(Xint) %*% epsilon  
  
  return(-as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  G <- gradientO(coef = coef, Y = Y, X = X)
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  hessian_sigma2 <- + n/(2*sigma2^2) - 2/(2*sigma2^3) * t(epsilon) %*% epsilon
  hessian_sigma2FDbeta <- - 1/(sigma2^2) * t(Xint) %*% epsilon
  hessian_beta <- - 1/(sigma2) * t(Xint) %*% Xint
  
  H <- cbind( rbind(hessian_beta,t(hessian_sigma2FDbeta)),
              rbind(hessian_sigma2FDbeta,hessian_sigma2)
  )
  attr(H,"grad") <- G
  
                        
  return(-H)
}



##### Gold standard
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seqPark_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seqParkNorm_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))

#### Internal derivatives
plvm.EPSODE0 <- estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = TRUE, trace = TRUE,
                        control = list(constrain = FALSE))

plvm.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2, 
                        fixSigma = FALSE, trace = TRUE, stepLambda1 = 10,
                        control = list(constrain = FALSE))

## check
plvm.test <- estimate(plvm.model, data = df.data, lambda1 = plvm.EPSODE$opt$message[3,"lambda1"])


#### External derivatives
external.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = TRUE, trace = TRUE,
                            objective = objectiveO, gradient = gradientO, hessian = hessianO, 
                            control = list(constrain = FALSE))

# graphical model: derivatives dSigma^-1 = Omega^-1 dOmega Omega^-1
