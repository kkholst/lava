rm(list = ls())

library(penalized)
library(lava)
library(testthat)
library(deSolve)
path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #### set the local path to the R files
source(file.path(path.lava,"tests","FCT.R"))
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

context("Reg-ElasticNet")

#### simulation ###
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

### models ###
lambda2 <- 50
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", lambda2 = lambda2, trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))


StepPenalize <- stepfun(x = rev(seq_lambda),
                        y = c(0,rev(unlist(lapply(penalized.PathL1, function(x){sum(x@penalized==0)} )))))
# curve(StepPenalize,0,2000)
#indexJump <- sapply(0:5, function(x){which(x ==  unlist(lapply(penalized.PathL1, function(x){sum(x@penalized==0)} )))[1]})


#### lasso path ###

#### LARS

test_that("LVM(LARS) vs penalize with ridge regression", {
  elvm.PathL1_LARS <- estimate(plvm.model,  data = df.data, lambda2 = lambda2, 
                               regularizationPath = 1)
  pathTempo <- penPath(elvm.PathL1_LARS)[,"lambda1.abs"]
  
  z = apply(outer(pathTempo,seq_lambda,'-'), 1, function(x){which.min(abs(x))})
  expect_equal(pathTempo, seq_lambda[z], tolerance=1e-2, scale= NULL)  
})

#### EPSODE
test_that("LVM(EPSODE-forward) vs penalize with ridge regression", {
  elvm.PathL1_EPSODE.f <- estimate(plvm.model,  data = df.data, regularizationPath = 2, lambda2 = lambda2, trace = TRUE)
  pathTempo <- penPath(elvm.PathL1_EPSODE.f)[,"lambda1.abs"]
  
  z = apply(outer(pathTempo,seq_lambda,'-'), 1, function(x){which.min(abs(x))})
  expect_equal(pathTempo, seq_lambda[z], tolerance=1e-2, scale= NULL)  
})


test_that("LVM(EPSODE-backward) vs penalize with ridge regression", {
  elvm.PathL1_EPSODE.b <- estimate(plvm.model,  data = df.data, increasing = FALSE, lambda2 = lambda2,
                                 regularizationPath = 2, trace = TRUE)
  pathTempo <- penPath(elvm.PathL1_EPSODE.b)[-nrow(penPath(elvm.PathL1_EPSODE.b)),"lambda1.abs"]
  
  z = apply(outer(pathTempo,seq_lambda,'-'), 1, function(x){which.min(abs(x))})
  expect_equal(pathTempo, seq_lambda[z], tolerance=1e-2, scale= NULL)  
})

#### check fix lambda ####
beta.free <- NULL
beta.fixed <- NULL

for(iter_l in 1:length(seq_lambda)){
  cat("*")
  eplvm.fit_tempo1 <- estimate(plvm.model,  data = df.data, fixSigma = FALSE, 
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2, 
                               lambda2 = lambda2/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                               control = list(constrain = TRUE))
  
   # normal model
  test_that("LVM vs pLVM with lasso", {
    expect_equal(object=unname(validLVM(eplvm.fit_tempo1, penalized.PathL1[[iter_l]])),
                 expected=rep(0,length(coef(eplvm.fit_tempo1))),
                 tolerance=0.01,scale=1)    
  })
  
  eplvm.fit_tempo2 <- estimate(plvm.model,  data = df.data, fixSigma = TRUE, 
                               lambda1 = penalized.PathL1[[iter_l]]@lambda1, 
                               lambda2 = lambda2, 
                               control = list(trace = FALSE))
  
  # fixed sigma
  test_that("LVM vs pLVM with lasso", {
    expect_equal(object=unname(validLVM(eplvm.fit_tempo2, penalized.PathL1[[iter_l]])),
                 expected=rep(0,length(coef(eplvm.fit_tempo2))),
                 tolerance=0.001,scale=1)    
  })
  
}



