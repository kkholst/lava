rm(list = ls())

library(penalized)
library(lava)
library(testthat)
library(deSolve)
path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #### set the local path to the R files
source(file.path(path.lava,"tests","FCT.R"))
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

context("Reg-lasso")

#### > standard regression ####

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
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

eplvm.model <- estimate(plvm.model, df.data, lambda1 = 40, lambda2 = 0)

StepPenalize <- stepfun(x = rev(seq_lambda),
                          y = c(0,rev(unlist(lapply(penalized.PathL1, function(x){sum(x@penalized==0)} )))))
# curve(StepPenalize,0,2000)
#indexJump <- sapply(0:5, function(x){which(x ==  unlist(lapply(penalized.PathL1, function(x){sum(x@penalized==0)} )))[1]})

#### no penalty ####
test_that("LVM vs pLVM with lasso - lambda=0", {
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, lambda2 = 0)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model),tolerance=0.001,scale=1)    
})

#### lasso path ####

#### LARS

test_that("LVM(LARS) vs penalize with lasso", {
  elvm.PathL1_LARS <- estimate(plvm.model,  data = df.data, 
                               regularizationPath = 1)
  indexJump_LARS <- sapply(0:5, function(x){which(x == rowSums(getPath(elvm.PathL1_LARS)[,names(coef(elvm.PathL1_LARS))]==0))[1]})
  
  StepLARS <- stepfun(x = getPath(elvm.PathL1_LARS)$lambda1.abs,
                          y = c(0,rowSums(getPath(elvm.PathL1_LARS)[,names(coef(elvm.PathL1_LARS))]==0))
  )
  # curve(StepLARS,0,2000)
  
  z = apply(outer(getPath(elvm.PathL1_LARS)[,"lambda1.abs"],seq_lambda,'-'), 2, function(x){min(abs(x))})
  expect_equal(z, rep(0,length(z)),tolerance=0.01,scale=1)    
})

#### EPSODE
test_that("LVM(EPSODE-forward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, regularizationPath = 2, control = list(trace = TRUE))
  
  z = apply(outer(getPath(elvm.PathL1_EPSODE)[,"lambda1.abs"],seq_lambda,'-'), 1, function(x){min(abs(x))})
  expect_equal(z, rep(0,length(z)),tolerance=0.01,scale=1)  
})

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
                                 regularizationPath = 2, control = list(trace = TRUE))
  
  z = apply(outer(getPath(elvm.PathL1_EPSODE)[-nrow(getPath(elvm.PathL1_EPSODE)),"lambda1.abs"],seq_lambda,'-'), 1, function(x){min(abs(x))})
  expect_equal(z, rep(0,length(z)),tolerance=0.01,scale=1)    
})

# lambda1.abs   lambda1 lambda2.abs lambda2             Y       Y~X1         Y~X2      Y~X3      Y~X4      Y~X5       Y,Y
# 7    0.000000   0.00000           0       0 -1.536800e-17 -0.0137584 -0.015866945 0.2384093 0.4476321 0.7105157 0.2110710
# 6    5.673622  26.79131           0       0 -1.538128e-17  0.0000000 -0.002845955 0.2265917 0.4357334 0.6993291 0.2117710
# 5    6.985398  32.94444           0       0 -1.536974e-17  0.0000000  0.000000000 0.2239570 0.4330683 0.6970158 0.2120358
# 4  129.031346 322.68978           0       0 -1.589213e-17  0.0000000  0.000000000 0.0000000 0.1998120 0.4637785 0.3998619
# 3  229.350069 334.18446           0       0 -1.396121e-17  0.0000000  0.000000000 0.0000000 0.0000000 0.2639665 0.6862978
# 2  361.069331 361.79292           0       0 -1.026874e-17  0.0000000  0.000000000 0.0000000 0.0000000 0.0000000 0.9980000
# 1  397.176622 397.97257           0       0 -1.026874e-17  0.0000000  0.000000000 0.0000000 0.0000000 0.0000000 0.9980000

#### check fix lambda ####
beta.free <- NULL
beta.fixed <- NULL

for(iter_l in 1:length(seq_lambda)){
  cat("*")
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                               control = list(constrain = TRUE, trace = FALSE))
  
  # normal model
  test_that("LVM vs pLVM with lasso", {
    expect_equal(object=unname(validLVM(eplvm.fit_tempo1, penalized.PathL1[[iter_l]])),
                 expected=rep(0,length(coef(eplvm.fit_tempo1))),
                 tolerance=0.001,scale=1)    
  })
 
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1)
  
  # fixed sigma
  test_that("LVM vs pLVM with lasso", {
    expect_equal(object=unname(validLVM(eplvm.fit_tempo2, penalized.PathL1[[iter_l]])),
                 expected=rep(0,length(coef(eplvm.fit_tempo1))),
                 tolerance=0.001,scale=1)    
  })
 
}

#### automatic estimation of sigma
source(file.path(path.lava,"tests","FCT.R"))
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})


elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, 
                               regularizationPath = 2, 
                               control = list(trace = TRUE))
  
  
#### > high dimensional ####

#### simulation ####
set.seed(10)
n <- 30
n.region <- 50
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

penalized.PathL1 <- penalized(Y ~  ., data = df.data, lambda1 = 1, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))


#### check fix lambda ####
beta.free <- NULL
beta.fixed <- NULL

# iter_l <- 5
# eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
#                              lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
#                              control = list(constrain = TRUE, trace = TRUE))
# 
# eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
#                              lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
#                              objective = objectiveO, gradient = gradientO, hessian = hessianO, 
#                              control = list(constrain = TRUE, trace = TRUE, abs.tol = 1e-10, rel.tol = 1e-9))


for(iter_l in 1:length(seq_lambda)){
  cat("*")
 
#   eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
#                                lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
#                                objective = objectiveLV, gradient = gradientLV, hessian = hessianLV, 
#                                control = list(constrain = TRUE, trace = TRUE, abs.tol = 1e-10, rel.tol = 1e-9))
#   
#   # normal model
#   test_that("penalized vs pLVM with lasso (high dimensional - sigmaFree)", {
#         expect_equal(object=unname(eplvm.fit_tempo1$par),
#                      expected=unname(coef2.penalized( penalized.PathL1[[iter_l]])),
#                      tolerance=1e-5)    
# #     expect_equal(object=unname(validLVM(eplvm.fit_tempo1, penalized.PathL1[[iter_l]])),
# #                  expected=rep(0,length(coef(eplvm.fit_tempo1))),
# #                  tolerance=0.01,scale=1)    
#   })
  
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                               objective = objectiveLV, gradient = gradientLV, hessian = hessianLV, 
                               control = list(constrain = FALSE, trace = FALSE, abs.tol = 1e-10, rel.tol = 1e-9))
  
  # fixed sigma
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFixed)", {
    expect_equal(object=unname(eplvm.fit_tempo2$par)[-length(eplvm.fit_tempo2$par)],
                 expected=unname(coef2.penalized( penalized.PathL1[[iter_l]]))[-length(eplvm.fit_tempo2$par)],
                 tolerance=1e-5)    
  })
}


