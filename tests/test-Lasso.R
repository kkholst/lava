rm(list = ls())

library(penalized)
library(lava)
library(testthat)
library(deSolve)
library(butils)
package.source("lava", Rcode = TRUE)
source(file.path(butils::dir.gitHub(),"lava","tests","FCT.R"))

context("Reg-lasso")

#### > standard regression ####
test.tolerance <- 1e-4
test.scale <- NULL

#### simulation ####
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

### models ####
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data, estimator = "gaussian") # undebug(lava:::estimate.lvm)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE, lambda2 = 0)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

#### no penalty ####
test_that("LVM vs pLVM with lasso - lambda=0", {
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, lambda2 = 0)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
  
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, lambda2 = 0, fixSigma = TRUE)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
})

#### lasso path ####

#### LARS

test_that("LVM(LARS) vs penalize with lasso", {
  elvm.PathL1_LARS <- estimate(plvm.model,  data = df.data, 
                               regularizationPath = 1, lambda2 = 0)
  
  coef0 <- unlist(getPath(elvm.PathL1_LARS, getCoef = "coef0", getLambda = NULL))
  indexJump_LARS <- sapply(0:5, function(x){which(x == coef0)[1]})
  
  lambda1path <- getLambda(elvm.PathL1_LARS, lambda1 = TRUE, abs = TRUE)[,1]
  StepLARS <- stepfun(x = lambda1path, y = c(0,coef0))
  # curve(StepLARS,0,2000)
  
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)    
})

#### EPSODE
test_that("LVM(EPSODE-forward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = TRUE,
                                 regularizationPath = 2, lambda2 = 0,
                                 control = list(trace =TRUE))
  
  lambda1path <- getLambda(elvm.PathL1_EPSODE, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)    
})
# getPath(elvm.PathL1_EPSODE)

# lambda1.abs   lambda1 lambda2.abs lambda2             Y          Y~X3          Y~X4          Y~X5       Y,Y
# 1    0.000000   0.00000           0       0 -1.578256e-17  2.384096e-01  4.476322e-01  7.105163e-01 0.2110710
# 2    5.674446  26.79518           0       0 -1.579584e-17  2.265903e-01  4.357318e-01  6.993281e-01 0.2117711
# 3    6.986304  32.94868           0       0 -1.578430e-17  2.239554e-01  4.330664e-01  6.970147e-01 0.2120359
# 4  129.032092 322.69010           0       0 -1.630669e-17 -1.289890e-06  1.998105e-01  4.637776e-01 0.3998638
# 5  229.350602 334.18454           0       0 -1.437577e-17  0.000000e+00 -1.126323e-06  2.639660e-01 0.6862993
# 6  361.070618 361.79319           0       0 -1.068329e-17  0.000000e+00  0.000000e+00 -1.942671e-06 0.9980028

# Regularization path: 
#   lambda1.abs   lambda1 lambda2.abs lambda2             Y          Y~X1          Y~X2          Y~X3          Y~X4          Y~X5       Y,Y
# 1    0.000000   0.00000           0       0 -1.578256e-17 -1.375840e-02 -1.586743e-02  2.384096e-01  4.476322e-01  7.105163e-01 0.2110710
# 2    5.673998  26.79308           0       0 -1.579584e-17  9.116710e-07 -2.845574e-03  2.265912e-01  4.357327e-01  6.993290e-01 0.2117710
# 3    6.986221  32.94829           0       0 -1.578430e-17  0.000000e+00  1.349033e-06  2.239556e-01  4.330666e-01  6.970149e-01 0.2120359
# 4  129.031622 322.69011           0       0 -1.630669e-17  2.868722e-18  1.217280e-18 -3.795390e-07  1.998114e-01  4.637786e-01 0.3998623
# 5  229.350637 334.18455           0       0 -1.437575e-17  7.389134e-20 -1.128473e-19  7.909855e-20 -1.195496e-06  2.639660e-01 0.6862993
# 6  361.070018 361.79324           0       0 -1.068329e-17  7.349418e-19 -3.339281e-19  9.422789e-19  8.002648e-21 -6.970495e-07 0.9980010

test_that("LVM(EPSODE-backward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = FALSE,
                                 regularizationPath = 2, lambda2 = 0,
                                 control = list(trace = TRUE))
  
  lambda1path <- getLambda(elvm.PathL1_EPSODE, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)    
})
# lambda1.abs   lambda1 lambda2.abs lambda2             Y       Y~X1         Y~X2      Y~X3      Y~X4      Y~X5       Y,Y
# 7    0.000000   0.00000           0       0 -1.536800e-17 -0.0137584 -0.015866945 0.2384093 0.4476321 0.7105157 0.2110710
# 6    5.673622  26.79131           0       0 -1.538128e-17  0.0000000 -0.002845955 0.2265917 0.4357334 0.6993291 0.2117710
# 5    6.985398  32.94444           0       0 -1.536974e-17  0.0000000  0.000000000 0.2239570 0.4330683 0.6970158 0.2120358
# 4  129.031346 322.68978           0       0 -1.589213e-17  0.0000000  0.000000000 0.0000000 0.1998120 0.4637785 0.3998619
# 3  229.350069 334.18446           0       0 -1.396121e-17  0.0000000  0.000000000 0.0000000 0.0000000 0.2639665 0.6862978
# 2  361.069331 361.79292           0       0 -1.026874e-17  0.0000000  0.000000000 0.0000000 0.0000000 0.0000000 0.9980000
# 1  397.176622 397.97257           0       0 -1.026874e-17  0.0000000  0.000000000 0.0000000 0.0000000 0.0000000 0.9980000


#### check fix lambda - at the breakpoints ####
# eplvm.model <- estimate(plvm.model, df.data, lambda1 = 10, lambda2 = 0)
for(iter_l in 1:length(seq_lambda)){
  cat("*")
  
  # normal model
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               lambda2 = 0, control = list(trace = -1))
  
  test_that("LVM vs pLVM with lasso", {
    expect_equal(object = coef(eplvm.fit_tempo1),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)    
  })
  
  # with constrains
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               lambda2 = 0, control = list(constrain = TRUE, trace = -1))
  
  test_that("LVM vs pLVM with lasso (constrain)", {
    expect_equal(object = coef(eplvm.fit_tempo2),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)  
  })
  
  # fixed sigma
  eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1, 
                               lambda2 = 0, control = list(trace = -1))
  
  test_that("LVM vs pLVM with lasso (fix sigma)", {
    expect_equal(object = coef(eplvm.fit_tempo3),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)  
  })
  
}

#### check fix lambda - between breakpoints ####
seq2_lambda <- igraph::running.mean(seq_lambda, binwidth = 2)

for(iter_l in 1:length(seq2_lambda)){
  cat("*")
  
  penalized.L1 <- penalized(Y ~  ., data = df.data, lambda1 = seq2_lambda[iter_l], trace = FALSE, lambda2 = 0)
  
  # constrain
  eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = seq2_lambda[iter_l], 
                               lambda2 = 0, control = list(trace = -1))
  
  test_that("LVM vs pLVM with lasso (fix sigma - between knots)", {
    expect_equal(object = coef(eplvm.fit_tempo3),
                 expected = coef2.penalized(penalized.L1),
                 tolerance = test.tolerance, scale=1)  
  })
  
  # no constrain - can find unexpected solution
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data,
                              lambda1 = seq2_lambda[iter_l]/penalized.L1@nuisance$sigma2, 
                              lambda2 = 0, control = list(trace = -1))
  
  test_that("LVM vs pLVM with lasso (fix sigma - between knots)", {
    penalized.L1bis <- penalized(Y ~  ., data = df.data, lambda1 = eplvm.fit_tempo1$penalty$lambda1.abs, trace = FALSE, lambda2 = 0)
    expect_equal(object = coef(eplvm.fit_tempo1),
                 expected = coef2.penalized(penalized.L1bis),
                 tolerance = test.tolerance, scale=1)
    
    
  })
}
cat("\n")

#### partial penalization ####
penalized.PathL1 <- penalized(Y~.,data = df.data, steps = "Park", trace = TRUE, lambda2 = 0, 
                              unpenalized = ~ X2 + X4,
                              penalized = Y ~ X1 + X3 + X5)
# penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE, lambda2 = 0)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

plvm.model2 <- penalize(lvm.model, value = Y ~ X1 + X3 + X5) # penalize(lvm.model)

test_that("LVM(LARS) vs penalize with lasso", {
  elvm.PathL1_LARS <- estimate(plvm.model2,  data = df.data, fit = NULL,
                                 regularizationPath = 1, lambda2 = 0)
 
  lambda1path <- getLambda(elvm.PathL1_LARS, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance*100, scale=test.scale)    
})


test_that("LVM(EPSODE-forward) vs penalize with lasso", {
  
  system.time(
  elvm.PathL1_EPSODEf <- estimate(plvm.model2,  data = df.data, increasing = TRUE, fit = NULL, estimator = "penalized",
                                 regularizationPath = 2, lambda2 = 0,
                                  control = list(trace =TRUE))
  )
  lambda1path <- getLambda(elvm.PathL1_EPSODEf, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)  
})

test_that("LVM(EPSODE-backward) vs penalize with lasso", {
  elvm.PathL1_EPSODEb <- estimate(plvm.model2,  data = df.data, increasing = FALSE, fit = NULL,
                                 regularizationPath = 2, lambda2 = 0,
                                 control = list(trace =TRUE))
  
  lambda1path <- getLambda(elvm.PathL1_EPSODEb, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda[seq_lambda>0],'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda[seq_lambda>0], tolerance=test.tolerance, scale=test.scale)    
})

#### > high dimensional ####

#### simulation ####
set.seed(10)
n <- 20
n.region <- 25
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:n.region), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( rbinom(n.region, size = 1, prob = 0.3)*rnorm(n.region, 3,1) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

### models ####
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, lambda1 = 1, steps = "Park", lambda2 = 0,trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
# 14.582073 14.580614 12.173491  8.819824  6.379138  3.257475  2.841114  2.619117  1.539485  1.440283  1.440139  1.376496  1.376358  1.000000

#### check fix lambda ####
for(iter_l in 1:length(seq_lambda)){
  cat("*")
  
   # fixed sigma
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE, lambda2 = 0,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                               control = list(constrain = TRUE, trace = -1))
  
  # print(range(coef(eplvm.fit_tempo2)-coef2.penalized( penalized.PathL1[[iter_l]])))
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFixed)", {
    expect_equal(object=coef(eplvm.fit_tempo2),
                 expected=coef2.penalized( penalized.PathL1[[iter_l]]),
                 tolerance=1e-3)
  })
  
  # normal model - can find unexpected solution
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE, lambda2 = 0,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                               control = list(constrain = TRUE, trace = -1, abs.tol = 1e-10, rel.tol = 1e-9))
  
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFree)", {
    penalized.L1bis <- penalized(Y ~  ., data = df.data, lambda1 = eplvm.fit_tempo1$penalty$lambda1.abs, trace = FALSE, lambda2 = 0)
    expect_equal(object=coef(eplvm.fit_tempo1),
                 expected=coef2.penalized( penalized.L1bis ),
                 tolerance=1e-2)
  })
}

#### regularization path ####
test_that("LVM(EPSODE-backward) vs penalize with lasso (high dimensional)", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = FALSE,
                                 regularizationPath = 2, lambda2 = 0, stopLambda = min(seq_lambda)/1.05,
                                 control = list(trace = TRUE))
  lambda1path <- getLambda(elvm.PathL1_EPSODE, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda[-length(seq_lambda)],'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda[-length(seq_lambda)], tolerance=10*test.tolerance, scale=test.scale)    
  
  for(iter_l in 1:length(indexLambda)){
    index_l <- which.min(abs(elvm.PathL1_EPSODE$regularizationPath$lambda1.abs-seq_lambda[iter_l]))
    expect_equal(object=unlist(elvm.PathL1_EPSODE$regularizationPath[index_l,-(1:5)]),
                 expected=coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance=0.1, scale=test.scale)  
  }
  plot(elvm.PathL1_EPSODE, type = "path")
  
})



#### > factor analysis ####

#### Simulations ####
set.seed(10)
n <- 500
formula.lvm1 <- as.formula( paste0("Y1~eta+", paste0("X",1, collapse = "+") ) )
formula.lvm2 <- as.formula( paste0("Y2~eta+", paste0("X",2, collapse = "+") ) )
formula.lvm3 <- as.formula( paste0("Y3~eta+", paste0("X",3, collapse = "+") ) )
formula.lvm4 <- as.formula( paste0("Y4~eta+", paste0("X",4:5, collapse = "+") ) )
formula.lvm5 <- Y5~1

lvm.modelSim <- lvm(list(formula.lvm1,
                         formula.lvm2,
                         formula.lvm3,
                         formula.lvm4,
                         formula.lvm5))
distribution(lvm.modelSim,~eta) <- normal.lvm(sd = 2)
latent(lvm.modelSim) <- ~eta
df.data <- sim(lvm.modelSim,n)
df.data <- df.data[,names(df.data) != "eta"]

formula.All_lvm <- sapply(paste0("Y",1:5,"~eta+", paste0("X",1:5, collapse = "+") ), as.formula)

#### models ####
lvm.model <- lvm(formula.All_lvm)
latent(lvm.model) <- "eta"
plvm.model <- penalize(lvm.model)
plvm.modelSim <- penalize(lvm.modelSim)

#### no penalty ####
test_that("LVM vs pLVM (lambda = 0)", {
  resLVM <- estimate(lvm.modelSim, data = scale(df.data))
  
  resPLVM <- estimate(plvm.modelSim, data = df.data, lambda1 = 0)
  expect_equal(coef(resLVM), coef(resPLVM), tolerance=test.tolerance, scale=test.scale)
  
  resPLVM <- estimate(plvm.modelSim, data = df.data, lambda1 = 0, control = list(constrain = TRUE))
  expect_equal(coef(resLVM), coef(resPLVM), tolerance=test.tolerance, scale=test.scale)  
  
  resPLVM <- estimate(plvm.modelSim, data = df.data, lambda1 = 0, fixSigma = TRUE)
  expect_equal(coef(resLVM), coef(resPLVM), tolerance=test.tolerance, scale=test.scale)
})

#### A given sigma ####
lambda1 <- 67.24619
test_that("LVM vs pLVM (lambda > 0)", {
  resPLVM1 <- estimate(plvm.modelSim, data = df.data, lambda1 = lambda1, lambda2 = 0)
  resPLVM2 <- estimate(plvm.modelSim, data = df.data, lambda1 = lambda1, lambda2 = 0, control = list(constrain = TRUE))
  expect_equal(coef(resPLVM1), coef(resPLVM2), tolerance=test.tolerance, scale=test.scale)  
  
  norm <- sum(coef(resPLVM2)[paste(endogenous(resPLVM2),endogenous(resPLVM2),sep = ",")])
  resPLVM3 <- estimate(plvm.modelSim, data = df.data, lambda1 = lambda1*norm, lambda2 = 0, fixSigma = TRUE)
  expect_equal(coef(resPLVM2), coef(resPLVM3), tolerance=test.tolerance, scale=test.scale)
})

#### IN PROGRESS ####
test <- FALSE
if(test){
  
  #### regression
  # simul
  set.seed(10)
  n <- 500
  formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
  lvm.modelSim <- lvm()
  regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
  distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
  df.data <- sim(lvm.modelSim,n)
  df.data <- as.data.frame(scale(df.data))
  
  # models 
  lvm.model <- lvm(formula.lvm)
  elvm.model <- estimate(lvm.model, df.data)
  plvm.model <- penalize(lvm.model)
  
  penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE, lambda2 = 0)
  seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

  # derivatives
  lvGaussianO <- function(x){lvGaussian(x, df.data[[1]], as.matrix(df.data[,-1]))}
  scoreGaussianO <- function(x){scoreGaussian(x, df.data[[1]], as.matrix(df.data[,-1]))}
  hessianGaussianO <- function(x){hessianGaussian(x, df.data[[1]], as.matrix(df.data[,-1]))}
  
  # normal
  system.time(
    P1 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, estimator = "penalized",
                   regularizationPath = 2, lambda2 = 0,
                   control = list(trace =TRUE))
  )
  
  system.time(
    P2 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, estimator = "penalized",
                   regularizationPath = 2, lambda2 = 0,
                   objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO, 
                   control = list(trace =TRUE))
  )
  
  lambda1path <- getLambda(P1, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)  
  
  
  # check why penalized1 give different results compared to lvGaussian0
  P3 <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, fixSigma = FALSE,
                 regularizationPath = 2, lambda2 = 0, resolution_lambda1 = c(1e-3,1e-2), estimator = "penalized1",
                 # objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO,
                 control = list(trace =TRUE, constrain = FALSE))
  
  getPath(P1)
  P3$message[!is.na(P3$message[,"indexChange"]),]
  
  P3bis <- estimate(plvm.model,  data = df.data, increasing = TRUE, fit = NULL, fixSigma = FALSE,
                    regularizationPath = 2, lambda2 = 0, resolution_lambda1 = c(1e-3,1e-4),
                    objective = lvGaussianO, gradient = scoreGaussianO, hessian = hessianGaussianO,
                    control = list(trace =TRUE, constrain = FALSE))
  
  estimate(plvm.model, data = df.data, lambda1 = )
  # > getPath(P1)
  # lambda1.abs   lambda1 lambda2.abs lambda2             Y          Y~X1          Y~X2          Y~X3          Y~X4          Y~X5       Y,Y
  # 1    0.000000   0.00000           0       0 -1.578256e-17 -1.375840e-02 -1.586743e-02  2.384096e-01  4.476322e-01  7.105163e-01 0.2110710
  # 2    5.674446  26.79518           0       0 -1.579584e-17  1.998178e-06 -2.844546e-03  2.265903e-01  4.357318e-01  6.993281e-01 0.2117711
  # 3    6.986304  32.94868           0       0 -1.578430e-17  0.000000e+00  1.585480e-06  2.239554e-01  4.330664e-01  6.970147e-01 0.2120359
  # 4  129.032092 322.69010           0       0 -1.630669e-17  0.000000e+00  0.000000e+00 -1.289890e-06  1.998105e-01  4.637776e-01 0.3998638
  # 5  229.350602 334.18454           0       0 -1.437577e-17  0.000000e+00  0.000000e+00  0.000000e+00 -1.126323e-06  2.639660e-01 0.6862993
  # 6  361.070618 361.79319           0       0 -1.068329e-17  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -1.942671e-06 0.9980028
  
  # > P3$message[!is.na(P3$message[,"indexChange"]),]
  # lambda1.abs   lambda1 lambda2.abs lambda2 indexChange             Y         Y~X1          Y~X2          Y~X3        Y~X4          Y~X5       Y,Y
  # 8           NA  26.80705           0      NA           2 -1.548527e-17 3.736763e-08 -2.846401e-03  2.265920e-01  0.43573345  6.993297e-01 0.2117706
  # 14          NA  32.96480           0      NA           3 -1.532813e-17 0.000000e+00  1.234041e-08  2.239568e-01  0.43306787  6.970160e-01 0.2120354
  # 21          NA 323.06255           0      NA           4 -6.379410e-18 0.000000e+00  0.000000e+00 -1.008235e-06  0.19981071  4.637778e-01 0.3996279
  # 74          NA 415.28443           0      NA           5 -1.322999e-17 0.000000e+00  0.000000e+00  0.000000e+00 -0.01609266  2.478744e-01 0.5720816
  # 81          NA 471.64577           0      NA           6 -9.314723e-18 0.000000e+00  0.000000e+00  0.000000e+00  0.00000000 -9.947705e-07 0.7661640
  
  #### factor model
  set.seed(10)
  n <- 500
  formula.lvm1 <- as.formula( paste0("Y1~eta+", paste0("X",1, collapse = "+") ) )
  formula.lvm2 <- as.formula( paste0("Y2~eta+", paste0("X",2, collapse = "+") ) )
  formula.lvm3 <- as.formula( paste0("Y3~eta+", paste0("X",3, collapse = "+") ) )
  formula.lvm4 <- as.formula( paste0("Y4~eta+", paste0("X",4:5, collapse = "+") ) )
  formula.lvm5 <- Y5~1
  
  lvm.modelSim <- lvm(list(formula.lvm1,
                           formula.lvm2,
                           formula.lvm3,
                           formula.lvm4,
                           formula.lvm5))
  distribution(lvm.modelSim,~eta) <- normal.lvm(sd = 2)
  latent(lvm.modelSim) <- ~eta
  df.data <- sim(lvm.modelSim,n)
  df.data <- df.data[,names(df.data) != "eta"]
  
  formula.All_lvm <- sapply(paste0("Y",1:5,"~eta+", paste0("X",1:5, collapse = "+") ), as.formula)
  
  # models
  lvm.model <- lvm(formula.All_lvm)
  latent(lvm.model) <- "eta"
  plvm.model <- penalize(lvm.model)
  plvm.modelSim <- penalize(lvm.modelSim)
  
  #### lasso path ####
  test <- FALSE
  if(test){
    PathLVM1_f <- estimate(plvm.modelSim, data = df.data, increasing = TRUE,
                                    regularizationPath = 2, lambda2 = 0,
                                    control = list(trace = TRUE))
    path <- getPath(PathLVM1_f, names = c("lambda1.abs",paste(endogenous(plvm.modelSim),endogenous(plvm.modelSim),sep = ",")))
    path
    
    ##
    ## ISSUE: this is incorrect because sigma is not update along the path
    ##
    
    getPath(PathLVM1_f)$lambda1^2/getPath(elvm.PathL1_EPSODEf)$lambda1.abs
    
    elvm.PathL1_EPSODEb <- estimate(plvm.modelSim, data = df.data, increasing = FALSE,
                                    regularizationPath = 2, lambda2 = 0,
                                    control = list(trace = TRUE))
    getPath(elvm.PathL1_EPSODEf)
    getPath(elvm.PathL1_EPSODEb)
    
    plot(getPath(elvm.PathL1_EPSODEf)$lambda1.abs)
    
    PathLVM1_f.free <- estimate(plvm.modelSim, data = df.data, increasing = TRUE, fixSigma = FALSE,
                                regularizationPath = 2, lambda2 = 0,
                                control = list(trace = TRUE))
    plot(getPath(PathLVM1_f.free)$lambda1.abs)
    # plot(PathLVM1_f.free, type = "path")
    
    test_that("LVM EPSODE-forward init", {
      res <- estimate(plvm.modelSim, data = df.data, lambda1 = 0)
      expect_equal(unlist(getPath(PathLVM1_f, row = 1, getLambda = NULL)),
                   coef(res), tolerance=test.tolerance, scale=test.scale)
      
      
      lambda1.tempo <- getPath(PathLVM1_f, row = 3, getLambda = c("lambda1","lambda1.abs"), getCoef = NULL)
      
      resFixed <- estimate(plvm.modelSim, data = df.data, fixSigma = FALSE,
                      lambda1 = lambda1.tempo$lambda1)
      resFixed.abs <- estimate(plvm.modelSim, data = df.data, fixSigma = TRUE,
                           lambda1 = resFixed$penalty$lambda1.abs)
      
      coefPath <- unlist(getPath(PathLVM1_f, row = 3, getLambda = NULL))
      coefPath-coef(resFixed)
    })
    
    
    
    test_that("LVM EPSODE-forward == EPSODE-backward", {
      expect_equal(getPath(elvm.PathL1_EPSODEf),
                   getPath(elvm.PathL1_EPSODEb))
    })
    
    
    # getPath(elvm.PathL1_EPSODEf, getCoef = "penalized")
    # getPath(elvm.PathL1_EPSODEb, getCoef = "penalized")
    
    
    test_that("LVM(EPSODE-forward) vs penalize with lasso", {
      elvm.PathL1_EPSODE <- estimate(plvm.modelSim, data = df.data, increasing = TRUE,
                                     regularizationPath = 2, lambda2 = 0, stopParam = 5,
                                     control = list(trace = TRUE))
      regPath <- getPath(elvm.PathL1_EPSODE, getLambda = c("lambda1.abs","lambda1"), rm.duplicated = TRUE)
      
      for(iterLambda in 1:length(regPath$lambda1.abs)){
        elvm.tempo <- estimate(plvm.modelSim, data = df.data,
                               lambda1 = regPath[iterLambda,"lambda1"], lambda2 = 0) # regPath[iterLambda,"lambda1.abs"]/regPath[iterLambda,"Y1,Y1"]
        
        print(rbind(Proximal = coef(elvm.tempo)[elvm.PathL1_EPSODE$penalty$names.penaltyCoef], 
                    EPSODE = regPath[iterLambda,elvm.PathL1_EPSODE$penalty$names.penaltyCoef]))
        print(range(coef(elvm.tempo)-regPath[iterLambda,-(1:2)]))
        cat("\n")
      }
      
      
      
    })
  }
  
  #### add covariance ####
  lvm.extended <- extendModel(lvm.modelSim, type = "all")
  plvm.extended <- penalize(lvm.extended, covariance = TRUE, latent = TRUE)
  test_that("LVM vs pLVM (lambda = 0)", {
    resLVM <- estimate(lvm.extended, data = df.data)
    
    resPLVM <- estimate(plvm.extended, data = df.data, lambda1 = 0, lambda2 = 0)
    
    resPLVM <- estimate(plvm.extended, data = df.data, lambda1 = 10, lambda2 = 0, control = list(constrain = TRUE))
    
    resPLVM <- estimate(plvm.extended, data = df.data, lambda1 = 1e3, lambda2 = 0, fixSigma = TRUE, control = list(constrain = TRUE))
    # expect_equal(coef(resLVM), coef(resPLVM), tolerance=test.tolerance, scale=test.scale)
  })
}