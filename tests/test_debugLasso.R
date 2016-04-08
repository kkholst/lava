rm(list = ls())


#### load functions ####

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

#### 1- regression ####
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)

###
elvm1.Path <- estimate(plvm.model,  data = df.data,
                       regularizationPath = TRUE, 
                       fix.sigma = TRUE,
                       control = list(constrain = TRUE, iter.max = 5000, data = df.data))

penalized.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = df.data[,all.vars(formula.lvm)[-1]], 
                           step = "Park", trace = FALSE)
TruePath <- matrix(unlist(lapply(penalized.fit,
                          function(x){
                            c(x@lambda1,x@unpenalized,x@penalized,x@nuisance$sigma2)
                          })),
                   ncol = ncol(elvm1.Path$opt$message), byrow = TRUE)
TruePath
elvm1.Path$opt$message


elvm1.Path <- estimate(plvm.model,  data = df.data,
                       regularizationPath = TRUE, 
                       fix.sigma = TRUE, use.lavaDeriv = TRUE,
                       control = list(constrain = TRUE, iter.max = 5000, data = df.data))
elvm1.Path$opt$message


### debug
dataY <- as.matrix(df.data[,1,drop = FALSE])
dataX <- as.matrix(df.data[,-1,drop = FALSE])

objectiveLv <- function(coef, ...){
  
  XX <- cbind(1,dataX)
  YY <- dataY
  
  beta <- coef[-length(coef)]
  residuals <- YY-XX %*% cbind(beta)
  sigma2 <- coef[length(coef)]
  
  lv <- - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(residuals) %*% residuals
  # pen <-  - lambda1 * sum( abs(beta) ) - lambda2 / 2 * sum( beta^2 )
  
  return(drop(lv))
}

gradientLv <- function(coef, ...){
  
  XX <- cbind(1,dataX)
  YY <- dataY
  
  beta <- coef[-length(coef)]
  residuals <- YY-XX %*% cbind(beta)
  sigma2 <- coef[length(coef)]
  
  grad_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * crossprod(residuals, residuals)
  grad_beta <- + 1/(sigma2) * crossprod(XX, residuals)
  
  return(drop(c(grad_beta, grad_sigma2)))
}

hessianLv <- function(coef, ...){
  
  XX <- cbind(1,dataX)
  YY <- dataY
  
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  residuals <- YY-XX %*% cbind(beta)
  
  hess_sigma2 <- + n/(2*sigma2^2) - 1/(sigma2^3) * crossprod(residuals, residuals)
  hess_sigma2FDbeta <- - 1/(sigma2^2) * crossprod(XX, residuals)
  hess_beta <- - 1/(sigma2) * crossprod(XX, XX)
  
  return(cbind( rbind(hess_beta,t(hess_sigma2FDbeta)),
                rbind(hess_sigma2FDbeta,hess_sigma2)
  ))
}


elvm1.Path <- estimate(plvm.model,  data = df.data,
                       regularizationPath = TRUE, 
                       fix.sigma = TRUE,
                       control = list(constrain = TRUE, iter.max = 5000, browser = TRUE, data = df.data))

objective(start)
objectiveLv(start)

gradient(start)
gradientLv(start)

hessian(start) # , type = "hessian"
hessianLv(start)
