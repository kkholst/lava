rm(list = ls())
### issue with data normalisation
### issue with gradiant from lava
### issue with FISTA
### no need to fix sigma with a TRUE LVM model

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
formula.lvm <- Y~X1+X2+X3+X4
lvm.model <- lvm(list(formula.lvm))
df.data <- sim(lvm.model,n)
plvm.model <- penalize(lvm.model)


lvm.fit <- estimate(lvm.model,  data = df.data,
                      control = list(constrain = TRUE))

#### check fix lambda ####
penalized.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = df.data[,all.vars(formula.lvm)[-1]], 
                           steps = "Park")

coef(penalized.fit[[4]])

elvm1ISTA.fit <- estimate(plvm.model,  data = df.data, lambda1 = penalized.fit[[4]]@lambda1, fix.sigma = TRUE,
                           control = list(constrain = TRUE, iter.max = 5000))
coef(elvm1ISTA.fit)
elvm1ISTA.fit$opt$iterations

elvm1FISTA.fit <- estimate(plvm.model,  data = df.data, lambda1 = penalized.fit[[4]]@lambda1, proxGrad.method = "FISTA", fix.sigma = TRUE,
                           control = list(constrain = TRUE, iter.max = 5000))
coef(elvm1FISTA.fit)
elvm1FISTA.fit$opt$iterations

#### check regularization path ####
# resLassoPath <- LassoPath(dataX = dataX_norm, 
#                           dataY = df.data[,all.vars(formula.lvm)[1]], 
#                           iter_max = 100)

df.data_norm <- normData(df.data, formula = formula.lvm)$data

p1lvm.fit <- estimate(plvm.model,  data = df.data_norm, regularizationPath = TRUE, fix.sigma = TRUE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data_norm))
p1lvm.fit$opt$message

penalized.fit <- penalized(response = df.data_norm[,all.vars(formula.lvm)[1]], 
                           penalized = df.data_norm[,all.vars(formula.lvm)[-1]], 
                           steps = "Park")
unlist(lapply(penalized.fit, function(x){x@lambda1}))

#### incorrect when data are not normogonalized

p1lvm.fit <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fix.sigma = TRUE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data))
p1lvm.fit$opt$message

penalized.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = df.data[,all.vars(formula.lvm)[-1]], 
                           steps = "Park")
unlist(lapply(penalized.fit, function(x){x@lambda1}))


#### 2- LVM ####

set.seed(10)
n <- 500
lvm.model2 <- lvm()
regression(lvm.model2, Y1 ~ eta) <- 1
regression(lvm.model2, Y2 ~ eta + X1 + X2) <- c(1,rep(0.3,2))
regression(lvm.model2, Y3 ~ eta + X3 + X4 + X5) <- c(1,rep(0.5,3))
regression(lvm.model2, eta ~ X6 + X7 + X8 + X9 + X10) <- rep(1,5)
df.data2 <- sim(lvm.model2, 1e3)
df.data2 <- df.data2[, names(df.data2) != "eta"]

### correct specification
lvm.model2bis <- lvm()
regression(lvm.model2bis) <- Y1 ~ eta
regression(lvm.model2bis) <- Y2 ~ eta + X1 + X2
regression(lvm.model2bis) <- Y3 ~ eta + X3 + X4 + X5
regression(lvm.model2bis) <- eta ~ X6 + X7 + X8 + X9 + X10
latent(lvm.model2bis) <- ~eta

lvm.fit2 <- estimate(lvm.model2bis,  data = df.data2,
                     control = list(constrain = TRUE))
coef(lvm.fit2)

plvm.model2bis <- penalize(lvm.model2bis)

elvm2ISTA.fit <- estimate(plvm.model2bis,  data = df.data2, lambda1 = 0,
                          control = list(constrain = TRUE, iter.max = 5000, proxGrad.method = "ISTA", fix.sigma = FALSE))
range(coef(lvm.fit2) - coef(elvm2ISTA.fit))

elvm2FISTA.fit <- estimate(plvm.model2bis,  data = df.data2, lambda1 = 0,
                          control = list(constrain = TRUE, iter.max = 5000, proxGrad.method = "FISTA", fix.sigma = FALSE))
range(coef(lvm.fit2) - coef(elvm2FISTA.fit))

### incorrect specification
lvm.model2bis <- lvm()
regression(lvm.model2bis) <- Y1 ~ eta
regression(lvm.model2bis) <- Y2 ~ eta + X1 + X2 + X3 + X4 + X5
regression(lvm.model2bis) <- Y3 ~ eta + X1 + X2 + X3 + X4 + X5
regression(lvm.model2bis) <- eta ~ + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
latent(lvm.model2bis) <- ~eta

lvm.fit2 <- estimate(lvm.model2bis,  data = df.data2,
                     control = list(constrain = TRUE))
coef(lvm.fit2)

plvm.model2bis <- penalize(lvm.model2bis)

elvm2ISTA.fit <- estimate(plvm.model2bis,  data = df.data2, lambda1 = 40, proxGrad.method = "ISTA", fix.sigma = FALSE,
                          control = list(constrain = TRUE, iter.max = 5000))
# elvm2ISTA.fit <- estimate(plvm.model2bis,  data = df.data2, lambda1 = 40, proxGrad.method = "ISTA", fix.sigma = TRUE,
#                           control = list(constrain = TRUE, iter.max = 5000))

coef(elvm2ISTA.fit)
range(coef(lvm.fit2) - coef(elvm2ISTA.fit))

dt.ggplot <- data.table(coef = names(coef(lvm.fit2)), 
                        lvm = coef(lvm.fit2), plvm = coef(elvm2ISTA.fit))

require(ggplot2)
ggbase <- ggplot(  melt(dt.ggplot, id.vars = "coef", variable.name = "model"), aes(x = coef, y = value, group = model, col = model))
ggbase + geom_point()


# elvm2ISTA.fit <- estimate(plvm.model2bis,  data = df.data2, lambda1 = 40, proxGrad.method = "FISTA", fix.sigma = FALSE,
#                           control = list(constrain = TRUE, iter.max = 5000))
# 
# coef(elvm2FISTA.fit)

#### regularization path
elvm2ISTA.RP <- estimate(plvm.model2bis,  data = df.data2, regularizationPath = TRUE, proxGrad.method = "ISTA", fix.sigma = FALSE,
                         control = list(constrain = TRUE, iter.max = 1000))


# #### check objective
# lvm.test <- lvm(Y ~ 1)
# df.test <- sim(lvm.test, 1e3)
# 
# calcObj <- function(Y){
#   meanY <- mean(Y)
#   sdY <- sd(Y)
#   n <- length(Y)
#   
#   ## theorical
#   Omega <- diag(sdY^2, 1, 1)
#   xi <- mean(Y)
#   
#   ## empirical
#   Sigma <- diag(sdY^2, 1, 1)
#   mu <- mean(Y)
#   
#   TT <- Sigma + tcrossprod(mu-xi)
#   
#   sum(dnorm(Y, mean = meanY, sd =sdY, log = TRUE))
#   
#   -n/2*log(2*pi) - n/2 * log(det(Omega)) - n/2 * tr(TT %*% solve(Omega))
#   
# }
# calcObj(df.test$Y)
# 
# logLik(lm(Y ~ 1, data = df.test))
# 
# lvm.fitTest <- estimate(lvm.test,  data = df.test,
#                         control = list(constrain = TRUE))
# 
# lvm.fitTest$opt$objective
# 
# 
# lvm.simple <- lvm(Y1 ~ eta, Y2 ~ eta, Y3 ~ eta)
# latent(lvm.simple) <- ~ eta
# df.simple <- sim(lvm.simple,1e3)
# 
# fit.simple <- estimate(lvm.simple, df.simple[names(df.simple) != "eta"])
# coef(fit.simple)
