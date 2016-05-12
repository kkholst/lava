rm(list = ls())

#### set the path to the R files ####
path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"

#### loading ####
library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(deSolve)
library(lava)

vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

## log likelihood and derivative for a linear model
objectiveO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  lv <-  - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(epsilon) %*% epsilon
  
  return(as.numeric(lv))
}

gradientO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * t(epsilon) %*% epsilon    
  gradient_beta <- + 1/(sigma2) * t(Xint) %*% epsilon  
  
  return(as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
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
  
  return(H)
}

#### simulation ###
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

#### Orthogonalize data
# df.dataH <- df.data # df.data <- df.dataH
# resOrtho <- prepareDataPath.lvm(model = lvm.model, data = df.data, penalty = plvm.model$penalty)
# df.data[] <- resOrtho$data[]

#### models ####

## lvm
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

## regularization path - penalized package
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", lambda2 = 5, trace = TRUE)
seqPark_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seqParkNorm_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))

## regularization path - glmPath algorithm
plvm.glmPath <- estimate(plvm.model, data = df.data, regularizationPath = 1, fixSigma = TRUE)
plvm.glmPath

plvm.glmPath <- estimate(plvm.model, data = df.data, lambda1 = 6.286164, trace = TRUE, fixSigma = TRUE,
                         control = list(iter.max = 100, constrain = TRUE))


# 1 86.153347       0 -0.121409118  0.000000e+00  0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000 18.778273
# 2 86.326405       0 -0.121386816  0.000000e+00  0.000000e+00 0.000000e+00 0.000000e+00 0.0003032559 18.738755
# 3 78.974093       0 -0.036645741  0.000000e+00  0.000000e+00 0.000000e+00 4.628701e-09 1.1525888808 12.699449
# 4 77.484400       0  0.007523891  0.000000e+00  0.000000e+00 5.377920e-09 7.754719e-01 1.8912276217  7.826057
# 5  7.398226       0  0.052847809  0.000000e+00  0.000000e+00 9.050240e-01 1.861319e+00 2.9244174862  3.980388
# 6  7.397491       0  0.052848212  0.000000e+00 -6.790435e-06 9.050291e-01 1.861325e+00 2.9244221768  3.980386
# 7  6.286164       0  0.053473845 -1.692839e-08 -1.026122e-02 9.127400e-01 1.870240e+00 2.9317331449  3.976663
# 8  0.000000       0  0.056772636 -5.865776e-02 -7.109498e-02 9.577815e-01 1.921950e+00 2.9777820863  3.963549

#### EPSODE regularization path 

### from no penalization
plvm.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2, lavaDerivatives = TRUE,
                    control = list(constrain = FALSE, step_lambda1 = 10, 
                                   start = coef(estimate(lvm.model, data = df.data))))

plvm.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2)

plvm.EPSODE <- estimate(plvm.model, data = df.data, regularizationPath = 2, lavaDerivatives = FALSE, fixSigma = TRUE, correctionStep = FALSE,
                        control = list(constrain = FALSE, step_lambda1 = 10, 
                                       start = coef(estimate(lvm.model, data = df.data))))


plvm.EPSODE2 <- estimate(plvm.model, data = df.data, regularizationPath = 2,
                         control = list(constrain = FALSE, step_lambda1 = 10, correction.step = FALSE,
                                        start = coef(estimate(lvm.model, data = df.data))))

### from full penalization
plvm.EPSODE.inv <- estimate(plvm.model, data = df.data, regularizationPath = 2,
                            control = list(constrain = FALSE, step_lambda1 = -10, iter.max = 1000))

# plvm.EPSODE.inv2 <- estimate(plvm.model, data = df.data, regularizationPath = 2,
#                              control = list(constrain = FALSE, step_lambda1 = -10, correction.step = FALSE, iter.max = 1000))


