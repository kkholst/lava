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
lambda2 <- 0

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

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE, lambda2 = lambda2)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

#### no penalty ####
test_that("LVM vs pLVM with lasso - lambda=0", {
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, lambda2 = lambda2)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model), tolerance=test.tolerance, scale=test.scale)    
})

#### lasso path ####

#### LARS

test_that("LVM(LARS) vs penalize with lasso", {
  elvm.PathL1_LARS <- estimate(plvm.model,  data = df.data, 
                               regularizationPath = 1, lambda2 = lambda2)
  
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
                                 regularizationPath = 2, lambda2 = lambda2,
                                 control = list(trace =TRUE))
  
  lambda1path <- getLambda(elvm.PathL1_EPSODE, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance*100, scale=test.scale)    
})

# Estimation of the nuisance parameter   lambda1.abs lambda1 lambda2.abs lambda2 indexChange             Y         Y~X1          Y~X2          Y~X3          Y~X4          Y~X5 Y,Y
# 1    0.000000      NA           0      NA          NA -1.578256e-17 -1.37584e-02 -1.586743e-02  2.384096e-01  4.476322e-01  7.105163e-01   1
# 2    5.673998      NA           0      NA           2 -1.579584e-17  9.11671e-07 -2.845574e-03  2.265912e-01  4.357327e-01  6.993290e-01   1
# 3    6.986221      NA           0      NA           3 -1.578430e-17  0.00000e+00  1.349033e-06  2.239556e-01  4.330666e-01  6.970149e-01   1
# 4  129.031622      NA           0      NA           4 -1.630669e-17  0.00000e+00  0.000000e+00 -3.795390e-07  1.998114e-01  4.637786e-01   1
# 5  229.350637      NA           0      NA           5 -1.437575e-17  0.00000e+00  0.000000e+00  0.000000e+00 -1.195496e-06  2.639660e-01   1
# 6  361.070018      NA           0      NA           6 -1.068329e-17  0.00000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -6.970495e-07   1
# - done 

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
                                 regularizationPath = 2, lambda2 = lambda2,
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

#### check fix lambda ####
beta.free <- NULL
beta.fixed <- NULL

for(iter_l in 1:length(seq_lambda)){
  cat("*")
  
  # normal model
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               lambda2 = lambda2)
  
  test_that("LVM vs pLVM with lasso", {
    expect_equal(object = coef(eplvm.fit_tempo1),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)    
  })
 
  # with constrains
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               lambda2 = lambda2, control = list(constrain = TRUE))
  
  test_that("LVM vs pLVM with lasso (constrain)", {
    expect_equal(object = coef(eplvm.fit_tempo2),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)  
  })
  
  # fixed sigma
  eplvm.fit_tempo3 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1, 
                               lambda2 = lambda2)
  
  test_that("LVM vs pLVM with lasso (fix sigma)", {
    expect_equal(object = coef(eplvm.fit_tempo3),
                 expected = coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance = test.tolerance, scale=1)  
  })
 
}


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

penalized.PathL1 <- penalized(Y ~  ., data = df.data, lambda1 = 1, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
# 14.582073 14.580614 12.173491  8.819824  6.379138  3.257475  2.841114  2.619117  1.539485  1.440283  1.440139  1.376496  1.376358  1.000000

#### check fix lambda ####
beta.free <- NULL
beta.fixed <- NULL

iter_l <- 5
eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                             control = list(constrain = FALSE, trace = TRUE))

eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                             control = list(constrain = FALSE, trace = TRUE))
coef(eplvm.fit_tempo1)
penalized.PathL1[[iter_l]]@penalized
penalized.PathL1
eplvm.fit_tempo1
# 14.2225510078081 *  1.2100e+00
eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                             objective = objectiveO, gradient = gradientO, hessian = hessianO,
                             control = list(constrain = TRUE, trace = TRUE, abs.tol = 1e-10, rel.tol = 1e-9))


for(iter_l in 1:length(seq_lambda)){
  cat("*")
 
  # eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
  #                              lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
  #                              #objective = objectiveLV, gradient = gradientLV, hessian = hessianLV,
  #                              control = list(constrain = TRUE, trace = TRUE, abs.tol = 1e-10, rel.tol = 1e-9))
  # 
  # # normal model
  # test_that("penalized vs pLVM with lasso (high dimensional - sigmaFree)", {
  #       expect_equal(object=coef(eplvm.fit_tempo1),
  #                    expected=coef2.penalized( penalized.PathL1[[iter_l]]),
  #                    tolerance=1e-5)
  # })
  
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                               # objective = objectiveLV, gradient = gradientLV, hessian = hessianLV, 
                               control = list(constrain = FALSE, trace = FALSE, abs.tol = 1e-10, rel.tol = 1e-9))
  
  # fixed sigma
  test_that("penalized vs pLVM with lasso (high dimensional - sigmaFixed)", {
    expect_equal(object=coef(eplvm.fit_tempo2),
                 expected=coef2.penalized( penalized.PathL1[[iter_l]]),
                 tolerance=1e-1)  
  })
}

#### regularization path ####
test_that("LVM(EPSODE-backward) vs penalize with lasso (high dimensional)", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = FALSE,
                                 regularizationPath = 2, lambda2 = lambda2, stopLambda = min(seq_lambda)/1.05,
                                 control = list(trace = TRUE))
  lambda1path <- getLambda(elvm.PathL1_EPSODE, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda[-length(seq_lambda)],'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda[-length(seq_lambda)], tolerance=test.tolerance, scale=test.scale)    
  
  for(iter_l in 1:length(indexLambda)){
    index_l <- which.min(abs(elvm.PathL1_EPSODE$regularizationPath$lambda1.abs-seq_lambda[iter_l]))
    expect_equal(object=unlist(elvm.PathL1_EPSODE$regularizationPath[index_l,-(1:5)]),
                 expected=coef2.penalized(penalized.PathL1[[iter_l]]),
                 tolerance=1e-1)  
  }
  
})

plotPath(elvm.PathL1_EPSODE)


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
})

# test_that("LVM vs pLVM (lambda = 0, fixSigma)", {
#   resLVM <- estimate(lvm.modelSim, data = scale(df.data))
#   
#   resPLVM <- estimate(plvm.modelSim, data = df.data, lambda1 = 0, fixSigma = TRUE)
#   expect_equal(coef(resLVM), coef(resPLVM), tolerance=test.tolerance, scale=test.scale)
# })
## do not work for fixed sigma

#### lasso path ####
elvm.PathL1_EPSODEf <- estimate(plvm.modelSim, data = df.data, increasing = TRUE,
                               regularizationPath = 2, lambda2 = lambda2,
                               control = list(trace = TRUE))

elvm.PathL1_EPSODEb <- estimate(plvm.modelSim, data = df.data, increasing = FALSE,
                               regularizationPath = 2, lambda2 = lambda2,
                               control = list(trace = TRUE))

test_that("LVM EPSODE-forward == EPSODE-backward", {
  expect_equal(getPath(elvm.PathL1_EPSODEf),
               getPath(elvm.PathL1_EPSODEb))
})


# getPath(elvm.PathL1_EPSODEf, getCoef = "penalized")
# getPath(elvm.PathL1_EPSODEb, getCoef = "penalized")


test_that("LVM(EPSODE-forward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.modelSim, data = df.data, increasing = TRUE,
                                 regularizationPath = 2, lambda2 = lambda2, stopParam = 5,
                                 control = list(trace = TRUE))
  regPath <- getPath(elvm.PathL1_EPSODE, getLambda = c("lambda1.abs","lambda1"), rm.duplicated = TRUE)
  
  for(iterLambda in 1:length(regPath$lambda1.abs)){
    elvm.tempo <- estimate(plvm.modelSim, data = df.data,
                           lambda1 = regPath[iterLambda,"lambda1"], lambda2 = lambda2) # regPath[iterLambda,"lambda1.abs"]/regPath[iterLambda,"Y1,Y1"]
    
    print(rbind(Proximal = coef(elvm.tempo)[elvm.PathL1_EPSODE$penalty$names.penaltyCoef], 
                EPSODE = regPath[iterLambda,elvm.PathL1_EPSODE$penalty$names.penaltyCoef]))
    print(range(coef(elvm.tempo)-regPath[iterLambda,-(1:2)]))
    cat("\n")
  }
  
  
  elvm.tempo <- estimate(plvm.modelSim, data = df.data,
                         lambda1 =80.59349 , lambda2 = 0)
  
  lvm2 <- lvm(list(Y1~eta,Y2~eta,Y3~eta,Y4~eta,Y5~1))
  estimate(lvm2, data = df.data)
  plvm2 <- penalize(lvm2)
  estimate(lvm2, data = df.data, lambda1 = 0)
  elvm.PathL1_EPSODE2 <- estimate(plvm.modelSim, data = df.data, increasing = FALSE,
                                 regularizationPath = 2, lambda2 = lambda2, stopParam = 5,
                                 control = list(trace = TRUE))
  regPath2 <- getPath(elvm.PathL1_EPSODE2, getLambda = c("lambda1.abs","lambda1"), rm.duplicated = TRUE)
  
  
  regPath2[6,]
  coef(elvm.tempo)
  elvm.tempo <- estimate(plvm.modelSim, data = df.data, 
                         lambda1 = regPath2[iterLambda,"lambda1.abs"])
  elvm.tempo$penalty
  regPath2[iterLambda,"lambda1.abs"]/0.15095
    
  elvm.tempo <- estimate(plvm.modelSim, data = df.data, 
                         lambda1 = regPath2[iterLambda,"lambda1.abs"]/regPath2[iterLambda,"Y1,Y1"])
  
  
  getPath(elvm.PathL1_EPSODE, getLambda = "lambda1.abs", getCoef = "coefn0")
  
  
  elvm.tempo <- estimate(plvm.modelSim, data = df.data, increasing = TRUE,
                         lambda1 = regPath[iterLambda,"lambda1.abs"]/0.29359, lambda2 = lambda2) # regPath[iterLambda,"lambda1.abs"]/regPath[iterLambda,"Y1,Y1"]
  
  
  ggPath <- plotPath(elvm.PathL1_EPSODE)
  ggPath + coord_cartesian(xlim = c(154,175), ylim = c(0,1e-2))
  
})

getPath(elvm.PathL1_EPSODE, type = "penalized")
elvm.PathL1_EPSODE
plotPath(elvm.PathL1_EPSODE, add.point = TRUE, line.size = 1)  + coord_cartesian(xlim = c(154,175), ylim = c(0,1e-2))

getPath(elvm.PathL1_EPSODE, getLambda = "lambda1.abs", getCoef = "coefn0")



test_that("LVM(EPSODE-backward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.modelSim, data = df.data, increasing = FALSE,
                                 regularizationPath = 2, lambda2 = lambda2)
  
  lambda1path <- getLambda(elvm.PathL1_EPSODE, lambda1 = TRUE, abs = TRUE)[,1]
  indexLambda <- apply(outer(lambda1path,seq_lambda,'-'), 2, function(x){which.min(abs(x))})
  
  expect_equal(lambda1path[indexLambda], expected=seq_lambda, tolerance=test.tolerance, scale=test.scale)    
})



# rowSums(penPath(res.EPSODE_fixed, type = "coef") == 0)
penPath(res.EPSODE_free)[,1:2]
penPath(res.EPSODE_fixed)[,1:2]

iterLambda <- 30
lambda_tempo <- penPath(res.EPSODE_fixed, type = "lambda1", row = iterLambda)+0.1#res.EPSODE$opt$message[10,"lambda1.abs"],
system.time(
  res_free <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                       lambda1 = lambda_tempo,
                       control = list(constrain = TRUE, iter.max = 5000))
)
M <- cbind(coef(res_free), 
           # free = as.numeric(penPath(res.EPSODE_free, type = "coef", row = iterLambda)),
           fixed = as.numeric(penPath(res.EPSODE_fixed, type = "coef", row = iterLambda)))
# colSums(M==0)

range( coef(res_free) - penPath(res.EPSODE_fixed, type = "coef", row = iterLambda) )

lava:::print.lvmfit(res_free)

# lvm.modelC <- lvm.model
# coeff <- penPath(res.EPSODE, type = "coef")[15,]
# coefff <- coeff[grep("~", names(coeff), fixed = TRUE)]
# for(iter_p in 1:length(coefff)){
#   regression(lvm.modelC, as.formula(names(coefff)[iter_p])) <- coefff[iter_p]
# }
# coefff <- coeff[setdiff(1:length(coeff), grep("~|,", names(coeff)))]
# for(iter_p in 1:length(coefff)){
#   intercept(lvm.modelC, as.formula(paste0("~",names(coefff)[iter_p]))) <- coefff[iter_p]
# }
# 
# estimate(lvm.modelC, data = df.data)

lvm.model
?constrain
penPath(res.EPSODE, type = "coef")[1,]
res.EPSODE[,1:15]


res.EPSODE$opt$message[,c(1,2)]

estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = FALSE, stepLambda1 = 20, trace = TRUE,
         control = list(constrain = FALSE, iter.max = 5000))



coef(res_free)
cbind(res.EPSODE$opt$message[iterLambda,-(1:4),drop = TRUE], coef(res_free))

coef(res_free) - res.EPSODE$opt$message[10,names(coef(res_free))]
which(coef(res_free) == 0)
which( res.EPSODE$opt$message[10,names(coef(res_free))] == 0)

rbind(coef(res_free), res.EPSODE$opt$message[10,names(coef(res_free))])


res_free <- estimate(plvm.model,  data = df.data, fixSigma = FALSE,
                     lambda1 = 10,
                     control = list(constrain = TRUE, iter.max = 5000))
res0 <- estimate(lvm.model,  data = df.data)
quantile(coef(res0) - coef(res_free))

#### 3- More complex LVM
lvm.modelSim2 <- lvm(list(Y1 ~ eta1,
                          Y2 ~ eta1,
                          Y3 ~ eta1,
                          Y4 ~ eta1,
                          Y5 ~ eta2,
                          Y6 ~ eta2,
                          Y7 ~ eta2,
                          X1 ~ eta3,
                          X2 ~ eta3,
                          X3 ~ eta3,
                          eta3 ~ eta1 + eta2))
latent(lvm.modelSim2) <- ~ eta3 + eta2 + eta1

df.data <- sim(lvm.modelSim2,n)
df.data <- df.data[,names(df.data) %in% paste0("eta",1:3) == FALSE]

lvm.model2 <- lvm(list(Y1 ~ eta1 + eta4,
                       Y2 ~ eta1 + eta4,
                       Y3 ~ eta1 + eta4,
                       Y4 ~ eta1,
                       Y5 ~ eta2,
                       Y6 ~ eta2,
                       Y7 ~ eta2,
                       X1 ~ eta3,
                       X2 ~ eta3,
                       X3 ~ eta3,
                       eta3 ~ eta1 + eta2 + eta4))
latent(lvm.model2) <- ~ eta3 + eta2 + eta1 + eta4

res3 <- estimate(lvm.model2,  data = df.data)
coef(res3)

plvm.model2 <- penalize(lvm.model2, value = c(paste0("Y",1:3,"~eta4"),"eta4,eta4","Y5~eta2","Y6~eta2"))

res3 <- estimate(plvm.model2,  data = df.data,
                 lambda1 = 3,
                 control = list(constrain = FALSE, iter.max = 1000, trace = TRUE, step = NULL))
coef(res3)

res3 <- estimate(plvm.model2,  data = df.data,
                 lambda1 = 2,
                 control = list(constrain = FALSE, iter.max = 1000, trace = TRUE))
coef(res3)

##
