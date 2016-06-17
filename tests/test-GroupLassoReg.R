rm(list = ls())

library(penalized)
library(lava)
library(testthat)
library(deSolve)
path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #### set the local path to the R files
source(file.path(path.lava,"tests","FCT.R"))
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

context("Reg-GroupLasso")

#### data ####
library(gglasso)
data(bardet)
group1 <- rep(1, times = 5)#rep(1:5,each=1)
bardet$x <- bardet$x[,1:5]
df.bardet <- data.frame(scale(data.frame(bardet)))

#### models ####

formula_bardety <- as.formula(paste0("y ~ ", paste(names(df.bardet)[names(df.bardet)!="y"], collapse = "+")))
lvm.model_bardety <- lvm(formula_bardety)
lvm.fit_bardety <- estimate(lvm.model_bardety, data = df.bardet)

plvm.model_bardety <- penalize(lvm.model_bardety)
plvm.model_GL <- penalize(lvm.model_bardety) 
plvm.model_GL$penalty$group.penaltyCoef[] <- 1

gglasso.fit <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")
seq_lambda <- gglasso.fit$lambda

### no penalization

test_that("LVM vs pLVM with group lasso - lambda=0", {
  plvm.fit_GL <- estimate(plvm.model_GL, data = df.bardet, lambda1 = 0)
  expect_equal(object=coef(plvm.fit_GL),expected=coef(lvm.fit_bardety),tolerance=0.001,scale=NULL)    
})
  

#### group lasso ####
plvm.fit_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                        lambda1 = 68*seq_lambda[1] * nrow(df.bardet),
                        control = list(constrain = FALSE, iter.max = 1000))
coef(plvm.fit_GL)

plvm.fit_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                        lambda1 = 69*seq_lambda[1] * nrow(df.bardet),
                        control = list(constrain = FALSE, iter.max = 1000))
coef(plvm.fit_GL)

# Mcoef <- NULL
# for(iter in 1:length(seq_lambda)){
#   eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
#                              lambda1 = constrain = FALSE, lambda1 = lambda[iter] * nrow(df.bardet), 
#                              control = list(constrain = FALSE, iter.max = 1000,
#                                             start =if(iter>1){Mcoef[nrow(Mcoef),]}else{NULL})
#   )
#   Mcoef <- rbind(Mcoef,
#                  coef(eplvm.model_GL))
# }


#### OLD

mTEST <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls", lambda = 0)
coef(mTEST)

eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                           lambda1 = mTEST$lambda[1] * nrow(df.bardet))




iterLambda <- 100
eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet, fixSigma = TRUE,
                           lambda1 = m1$lambda[iterLambda] * nrow(df.bardet),
                           control = list(constrain = FALSE, iter.max = 1000))
coef(eplvm.model_GL)[grep("X",names(coef(eplvm.model_GL)), fixed = TRUE)] - m1$beta[,iterLambda]



#### agreement between lasso and grouped lasso when dealing with one parameter
plvm.model <- penalize(lvm.model_bardety, value = "bardety~X1")

eplvm.model <- estimate(plvm.model,  data = df.bardet,
                           lambda1 = 5,
                           control = list(constrain = TRUE, iter.max = 1000, trace = TRUE))
plvm.model

plvm.model_GL <- plvm.model
plvm.model_GL$penalty$group.penaltyCoef[] <- 1
eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet,
                           lambda1 = 5,
                           control = list(constrain = TRUE, iter.max = 1000, trace = TRUE))
coef(eplvm.model) - coef(eplvm.model_GL)

m1 <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")

m1$lambda * length(bardet$y) 

m1$b0

m1$beta
m1$npasses

set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)
lvm.test0 <- estimate(lvm.model,  data = df.data)
plvm.test0 <- estimate(plvm.model,  data = df.data)
coef(lvm.test0)

lvm.test1 <- estimate(lvm.model2,  data = df.bardet)
plvm.test1 <- estimate(penalize(lvm.model2),  data = df.bardet)
coef(lvm.test1)
names(df.bardet)

####
library(grplasso)
data(splice)

## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

## Fit a logistic model 
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(),
                       contrasts = contr, center = TRUE, standardize = TRUE)

####
library(grpreg)
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="grLasso")
plot(fit)









