rm(list = ls())


#### load functions ####

library(lava)
library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(deSolve)
library(regsem)
library(testthat)
library(butils)
package.source("lava", Rcode = TRUE)
source(file.path(butils::dir.gitHub(),"lava","tests","FCT.R"))

n <- 500

#### 1- Simulations ####
set.seed(10)
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

lvm.model <- lvm(formula.All_lvm)
latent(lvm.model) <- "eta"
plvm.model <- penalize(lvm.model)
plvm.modelSim <- penalize(lvm.modelSim)

#### 2- Estimations ####

test_that("LVM vs pLVM (lambda = 0)", {
  resLVM <- estimate(lvm.modelSim, data = scale(df.data))
  resPLVM <- estimate(plvm.modelSim, data = df.data, lambda1 = 0)
  
  expect_equal(coef(resLVM), coef(resPLVM), tolerance=1e-4, scale= NULL)  
})

# test_that("LVM vs pLVM (lambda = 0)", {
#   resLVM <- estimate(lvm.modelSim, data = scale(df.data))
#   resPLVM <- estimate(plvm.modelSim, data = df.data, lambda1 = 0, fixSigma = TRUE)
#   
#   expect_equal(coef(resLVM), coef(resPLVM), tolerance=1e-4, scale= NULL)  
# })

res.EPSODE_free <- estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = FALSE, trace = TRUE,
                       control = list(constrain = FALSE, iter.max = 5000))
res.EPSODE_fixed <- estimate(plvm.model, data = df.data, regularizationPath = 2, fixSigma = TRUE,  trace = TRUE,
                       control = list(constrain = FALSE, iter.max = 5000))

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


#### 99 - LAVA vs REGSEM ####

HS <- data.frame(scale(HolzingerSwineford1939[,7:15]))
mod <- '
f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
'
outt = cfa(mod,HS)
summary(outt)

#### REGSEM
fit0.regsem <- regsem(outt,
                      lambda=0,type="lasso", 
                      optMethod="nlminb", gradFun="ram")
summary(fit0.regsem)

fitLOW.regsem <- regsem(outt,
                        lambda=0.1,type="lasso", 
                        optMethod="nlminb", gradFun="ram")
summary(fitLOW.regsem)

fitHIGH.regsem <- regsem(outt,
                         lambda=0.35,type="lasso",
                         optMethod="nlminb", gradFun="ram")
summary(fitHIGH.regsem)


fitLOW.multi_optim <- multi_optim(outt,max.try=40,
                                  lambda=0.1,type="lasso", verbose = FALSE,
                                  optMethod="nlminb", gradFun="ram")
summary(fitLOW.multi_optim)

fitHIGH.multi_optim <- multi_optim(outt,max.try=40,
                                   lambda=0.35,type="lasso", verbose = FALSE,
                                   optMethod="nlminb", gradFun="ram")
summary(fitHIGH.multi_optim)

fitHIGH.multi_optimX <-  try(multi_optim(outt,max.try=40,
                                         lambda=0.35,type="lasso", verbose = FALSE,
                                         optMethod="optimx", gradFun="ram")
)

#### LAVA
linkLVM <- paste0(paste0("x",1:9),"~eta")
lvm.HS <- lvm(lapply(linkLVM, as.formula))
latent(lvm.HS) <- "eta"
plvm.HS <- penalize(lvm.HS,linkLVM[-1])

elvm.HS <- estimate(lvm.HS, data = HS)
eplvm0.HS <- estimate(plvm.HS, data = HS, lambda1 = 0, control = list(start = coef(elvm.HS)), fixSigma = FALSE )
eplvmLOW.HS <- estimate(plvm.HS, data = HS, lambda1 = 0.1*nrow(HS), control = list(start = coef(elvm.HS)), fixSigma = FALSE )
eplvmHIGH.HS <- estimate(plvm.HS, data = HS, lambda1 = 0.35*nrow(HS), control = list(start = coef(elvm.HS)), fixSigma = FALSE )

EPSODE.HS <- estimate(plvm.HS, data = HS, regularizationPath = 2, fixSigma = TRUE, stepLambda1 = 100, trace = TRUE)
EPSODE.HS$opt$message[,c(1:4,14:21)]

resTest <- estimate(plvm.HS,  data = HS, lambda1 = 39.39044)

# quantile(coef(elvm.HS) - coef(eplvm0.HS))
# quantile(coef(elvm.HS) - coef(eplvm0.fixedHS))


resTest <- estimate(plvm.HS,  data = HS, lambda1 = 100,  control = list(constrain = TRUE))
coef(resTest)

rbind(
  LAVA = coef(elvm.HS)[linkLVM[-1]],
  REGSEM = fit0.regsem$coefficients[1:(length(linkLVM)-1)]
)

rbind(
  LAVA = coef(eplvmLOW.HS)[linkLVM[-1]],
  REGSEM = fitLOW.regsem$coefficients[1:(length(linkLVM)-1)]
)

rbind(
  LAVA = coef(eplvmHIGH.HS)[linkLVM[-1]],
  REGSEM = fitHIGH.regsem$coefficients[1:(length(linkLVM)-1)]
)