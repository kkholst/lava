#%%%%%
#%
#% Comparison penalized LVM with penalized package for linear models
#%
#%%%%%

rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(lava)
library(deSolve)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #"C:/Users/hpl802/Downloads/lava-penalization/lava-penalization"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

coef2.penalized <- function(x, iter_lambda){
  
  if(is.list(x)){
    
    if(!missing(iter_lambda)){
      x <- x[[iter_lambda]]
    }else{
      res <- lapply(x, function(model){
        c(model@lambda1,
          model@lambda2,
          model@unpenalized, 
          model@penalized, 
          model@nuisance$sigma2)
      })
      
      Mres <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
      colnames(Mres) <- c("lambda1","lambda2",names(x[[1]]@unpenalized),names(x[[1]]@penalized),"sigma2")
      return(Mres)
    }
    
  } 
  
  coef <- c(x@unpenalized, 
            x@penalized, 
            x@nuisance$sigma2)
  return(coef)
}

validLVM <- function(x){
  library(penalized)
  
  lambda1 <- x$penalty$lambda1
  if(is.null(x$control$proxGrad$fixSigma)){lambda1 <- lambda1*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]}
  lambda2 <- x$penalty$lambda2
  
  resPenalized <- penalized(as.formula(paste0(endogenous(x),"~.")), 
                            lambda1 = lambda1, 
                            lambda2 = lambda2, 
                            data = x$data$model.frame, trace = FALSE)
  diffCoef <- coef(x) - coef2.penalized(resPenalized)
  
  return(diffCoef)
}

#### 0- simulation ####
set.seed(20)
n <- 100
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),0.25,0.5,0.75) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)


#### 1- only L1 ####

####
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

#### check fix lambda ####
beta.free <- NULL
beta.fixed <- NULL

for(iter_l in 1:length(seq_lambda)){
  cat(iter_l, " : ",penalized.PathL1[[iter_l]]@lambda1,"\n")
  
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                             control = list(constrain = TRUE, trace = FALSE))
  
  beta.free <- rbind(beta.free,
                     c(iter = elvm.fit_tempo$opt$iterations, coef(elvm.fit_tempo), 
                       range(validLVM(elvm.fit_tempo)))
  )
 
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1)
  
  beta.fixed <- rbind(beta.fixed,
                     c(iter = elvm.fit_tempo$opt$iterations, coef(elvm.fit_tempo), 
                       range(validLVM(elvm.fit_tempo)))
  )
}

#### check path ####
elvm.PathL1_fixed <- estimate(plvm.model,  data = df.data, 
                              regularizationPath = 1,
                              control = list(constrain = TRUE, iter.max = 5000))
elvm.PathL1_fixed$opt$message

# lambda1.abs   lambda1 lambda2          Y    Y~X1_penY    Y~X2_penY    Y~X3_penY    Y~X4_penY    Y~X5_penY      Y,Y
# 1     0.000000  0.000000       0 -0.1524526 5.336037e-02 7.180407e-02 1.975278e-01 7.702844e-01 9.961457e-01 4.871482
# 9     5.694868  1.164618       0 -0.1569026 5.335898e-06 2.225382e-02 1.362848e-01 6.886257e-01 9.185232e-01 4.889901
# 8     5.695438  1.164734       0 -0.1569030 0.000000e+00 2.224886e-02 1.362787e-01 6.886175e-01 9.185154e-01 4.889904
# 7     7.850693  1.600830       0 -0.1583222 0.000000e+00 1.525408e-15 1.148891e-01 6.585703e-01 8.871274e-01 4.904138
# 6    19.371842  3.855989       0 -0.1760131 0.000000e+00 0.000000e+00 1.931809e-05 5.020948e-01 7.187859e-01 5.023832
# 5    19.373779  3.856353       0 -0.1760161 0.000000e+00 0.000000e+00 0.000000e+00 5.020685e-01 7.187576e-01 5.023860
# 4    55.973662  9.694055       0 -0.2445822 0.000000e+00 0.000000e+00 0.000000e+00 7.677923e-05 2.251495e-01 5.774019
# 3    55.979260  9.694738       0 -0.2445927 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 2.250740e-01 5.774190
# 2    77.057332 12.687444       0 -0.2452850 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 8.226082e-05 6.073511
# 11   77.065039 12.688449       0 -0.2452852 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 6.073638
# elvm.PathL1_fixed <- estimate(plvm.model,  data = df.data, 
#                               regularizationPath = 1,
#                               control = list(constrain = FALSE, iter.max = 5000))
# elvm.PathL1_fixed$opt$message

elvm.EPSODE1_fixed <- estimate(plvm.model,  data = df.data, 
                              regularizationPath = 2, fixSigma = TRUE, stepLambda1 = 20, correctionStep = FALSE,
                              control = list(constrain = FALSE, iter.max = 5000))

elvm.EPSODE1_free <- estimate(plvm.model,  data = df.data, 
                             regularizationPath = 2, fixSigma = FALSE, stepLambda1 = 2,
                             control = list(constrain = FALSE, iter.max = 5000))

elvm.EPSODE1_fixed <- estimate(plvm.model,  data = df.data, 
                               regularizationPath = 2, fixSigma = TRUE, stepLambda1 = -20, correctionStep = FALSE,
                               control = list(constrain = FALSE, iter.max = 5000))

#### 2- L1 and L2 ####
penalized.PathL12 <- penalized(Y ~  ., data = df.data, steps = "Park", lambda2 = 100, trace = FALSE)
seq_lambda <- unlist(lapply(penalized.PathL12, function(x){x@lambda1}))

beta.free <- NULL ### not ok
beta.fixed <- NULL

for(iter_l in 1:length(seq_lambda)){
  cat(iter_l, " : ",penalized.PathL12[[iter_l]]@lambda1,"\n")
  
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
                             lambda1 = penalized.PathL12[[iter_l]]@lambda1/penalized.PathL12[[iter_l]]@nuisance$sigma2, 
                             lambda2 = penalized.PathL12[[iter_l]]@lambda2,
                             control = list(constrain = TRUE, trace = FALSE))
  
  beta.free <- rbind(beta.free,
                     c(iter = elvm.fit_tempo$opt$iterations, coef(elvm.fit_tempo), 
                       range(validLVM(elvm.fit_tempo)))
  )
  
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, fixSigma = TRUE,
                             lambda1 = penalized.PathL12[[iter_l]]@lambda1, lambda2 = penalized.PathL12[[iter_l]]@lambda2)
  
  beta.fixed <- rbind(beta.fixed,
                      c(iter = elvm.fit_tempo$opt$iterations, coef(elvm.fit_tempo), 
                        range(validLVM(elvm.fit_tempo)))
  )
}

#### check path ####
elvm.PathL12_fixed <- estimate(plvm.model,  data = df.data, 
                              regularizationPath = 1, lambda2 = 100, 
                              control = list(constrain = TRUE, iter.max = 5000))
elvm.PathL12_fixed$opt$message

elvm.PathL12_fixed <- estimate(plvm.model,  data = df.data, 
                              regularizationPath = 1, lambda2 = 100, 
                              control = list(constrain = FALSE, iter.max = 5000))
elvm.PathL12_fixed$opt$message

#### TO CHECK
# elvm.EPSODE12_fixed <- estimate(plvm.model,  data = df.data, 
#                               regularizationPath = 2, fixSigma = TRUE, stepLambda1 = 5.5, lambda2 = 100, 
#                               control = list(constrain = FALSE, iter.max = 5000))


#### 3- Group Lasso ####
# gglasso
# grplasso
# grpreg

####
library(gglasso)
data(bardet)
group1 <- rep(1, times = 5)#rep(1:5,each=1)
bardet$x <- bardet$x[,1:3]

# fit group lasso penalized least squares
m1 <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")
m1$b0
m1$beta
m1$npasses

df.bardet <- data.frame(bardet$y, bardet$x)
names(df.bardet) <- gsub(".","",names(df.bardet), fixed = TRUE)
names(df.bardet)[1] <- "Y"
formula_bardety <- as.formula(paste0("Y ~ ", paste(names(df.bardet)[-1], collapse = "+")))
lvm.model_bardety <- lvm(formula_bardety)
plvm.model_bardety <- penalize(lvm.model_bardety)


elvm.model_bardety <- estimate(lvm.model_bardety,  data = df.bardet)
eplvm.model_bardety <- estimate(plvm.model_bardety,  data = df.bardet, lambda1 = 0, fixSigma = TRUE)

eplvm.model_G1 <- estimate(plvm.model_G1,  data = df.bardet,
                           lambda1 = 0, fixSigma = TRUE, trace = TRUE, 
                           control = list(constrain = FALSE, iter.max = 1000))


set.seed(20)
n <- nrow(df.bardet)
lvm.modelSim <- lvm(formula_bardety)
regression(lvm.modelSim, formula_bardety) <- as.list( coef(lm(formula_bardety, df.bardet))[-1] )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = sd(lm(formula_bardety, df.bardet)$residuals))
df.data <- sim(lvm.modelSim,n)

elvm.model <- estimate(lvm.model_bardety,  data = df.data)
eplvm.model <- estimate(plvm.model_bardety,  data = df.data, lambda1 = 0, fixSigma = TRUE)



df.data[,3]  <- df.bardet[,3] 


coef(eplvm.model_G1) - coef(elvm.model_bardety) ### not really converged

eplvm.model_G1 <- estimate(plvm.model_G1,  data = df.bardet,
                           lambda1 = 100, fix.sigma = TRUE, 
                           control = list(constrain = TRUE, iter.max = 1000, trace = FALSE))
coef(eplvm.model_G1)


plvm.model_GL <- penalize(lvm.model_bardety) 
plvm.model_GL$penalty$group.penaltyCoef[] <- 1

eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet,
                          lambda1 = 0, fix.sigma = TRUE, 
                          control = list(constrain = TRUE, iter.max = 1000, start = coef(elvm.model_bardety), trace = TRUE))
coef(eplvm.model_GL) - coef(elvm.model_bardety) 

eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet,
                           lambda1 = 0, fix.sigma = TRUE, 
                           control = list(constrain = TRUE, iter.max = 10, trace = TRUE))
coef(eplvm.model_GL)### not really converged

eplvm.model_GL <- estimate(plvm.model_GL,  data = df.bardet,
                           lambda1 = m1$lambda[2] * nrow(df.bardet), fix.sigma = TRUE, 
                           control = list(constrain = FALSE, iter.max = 100, trace = TRUE))
coef(eplvm.model_GL)
m1$beta[,2]

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


