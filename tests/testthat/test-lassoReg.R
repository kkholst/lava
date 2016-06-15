rm(list = ls())

library(penalized)
library(lava)
library(testthat)
library(deSolve)

context("LVM-lassoRegression")


path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #### set the local path to the R files
source(file.path(path.lava,"tests","FCT.R"))
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

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
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))


StepPenalize <- stepfun(x = rev(seq_lambda),
                          y = c(0,rev(unlist(lapply(penalized.PathL1, function(x){sum(x@penalized==0)} )))))
# curve(StepPenalize,0,2000)
#indexJump <- sapply(0:5, function(x){which(x ==  unlist(lapply(penalized.PathL1, function(x){sum(x@penalized==0)} )))[1]})

#### no penalty ####
test_that("LVM vs pLVM with lasso", {
  eplvm.model <- estimate(plvm.model, df.data, lambda1 = 0, lambda2 = 0)
  expect_equal(object=coef(elvm.model),expected=coef(eplvm.model),tolerance=0.001,scale=1)    
})

#### lasso path ###

#### LARS

test_that("LVM(LARS) vs penalize with lasso", {
  elvm.PathL1_LARS <- estimate(plvm.model,  data = df.data, 
                               regularizationPath = 1)
  indexJump_LARS <- sapply(0:5, function(x){which(x == rowSums(penPath(elvm.PathL1_LARS)[,names(coef(elvm.PathL1_LARS))]==0))[1]})
  
  StepLARS <- stepfun(x = penPath(elvm.PathL1_LARS)$lambda1.abs,
                          y = c(0,rowSums(penPath(elvm.PathL1_LARS)[,names(coef(elvm.PathL1_LARS))]==0))
  )
  # curve(StepLARS,0,2000)
  
  z = apply(outer(penPath(elvm.PathL1_LARS)[,"lambda1.abs"],seq_lambda,'-'), 2, function(x){min(abs(x))})
  expect_equal(z, rep(0,length(z)),tolerance=0.01,scale=1)    
})

#### EPSODE
test_that("LVM(EPSODE-forward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, regularizationPath = 2, control = list(trace = TRUE))
  
  z = apply(outer(penPath(elvm.PathL1_EPSODE)[,"lambda1.abs"],seq_lambda,'-'), 1, function(x){min(abs(x))})
  expect_equal(z, rep(0,length(z)),tolerance=0.01,scale=1)  
})


test_that("LVM(EPSODE-backward) vs penalize with lasso", {
  elvm.PathL1_EPSODE <- estimate(plvm.model,  data = df.data, increasing = FALSE,
                                 regularizationPath = 2, trace = TRUE)
  
  z = apply(outer(penPath(elvm.PathL1_EPSODE)[-nrow(penPath(elvm.PathL1_EPSODE)),"lambda1.abs"],seq_lambda,'-'), 1, function(x){min(abs(x))})
  expect_equal(z, rep(0,length(z)),tolerance=0.01,scale=1)    
})

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




