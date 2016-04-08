#### note
# talk with klauss about gradient / hessian : hessian(beta, type = "expectation") vs hessian(beta, type = "hessian")
#


### issue with data normalisation
### issue with gradient from lava
### issue with FISTA


rm(list = ls())

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
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)

lvm.fit <- estimate(lvm.model,  data = df.data,
                      control = list(constrain = TRUE))

#### check fix lambda ####
penalized.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = df.data[,all.vars(formula.lvm)[-1]], 
                           lambda1 = 500)

elvm1ISTA.fit <- estimate(plvm.model,  data = df.data, 
                          lambda1 = penalized.fit@lambda1/penalized.fit@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 5000))

coef(elvm1ISTA.fit) - c(penalized.fit@unpenalized,penalized.fit@penalized,penalized.fit@nuisance$sigma2)

penalized.fit1000 <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = df.data[,all.vars(formula.lvm)[-1]], 
                           lambda1 = 1000)

elvm1.FIXED <- estimate(plvm.model,  data = df.data, lambda1 = penalized.fit1000@lambda1,
                        fix.sigma = TRUE, 
                        control = list(constrain = TRUE, iter.max = 5000))
coef(elvm1.FIXED) - c(penalized.fit1000@unpenalized, penalized.fit1000@penalized, penalized.fit1000@nuisance$sigma2)

elvm1.FREE <- estimate(plvm.model,  data = df.data, lambda1 = penalized.fit1000@lambda1/penalized.fit1000@nuisance$sigma2,
                       fix.sigma = FALSE,
                       control = list(constrain = TRUE, iter.max = 5000))
coef(elvm1.FREE) - c(penalized.fit1000@unpenalized, penalized.fit1000@penalized, penalized.fit1000@nuisance$sigma2)

####

penalizedPath.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                               penalized = df.data[,all.vars(formula.lvm)[-1]], 
                               steps = "Park")
iter_path <- 2

elvm1ISTA.fit <- estimate(plvm.model,  data = df.data, fix.sigma = FALSE,
                          lambda1 = penalizedPath.fit[[iter_path]]@lambda1,#/penalizedPath.fit[[iter_path]]@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 5000, cooling = NA))
elvm1ISTA.fit2 <- estimate(plvm.model,  data = df.data, fix.sigma = FALSE,
                          lambda1 = penalizedPath.fit[[iter_path]]@lambda1,#/penalizedPath.fit[[iter_path]]@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 5000, cooling = 0.999))
elvm1ISTA.fit$opt$iterations
elvm1ISTA.fit2$opt$iterations

coef(elvm1ISTA.fit) - c(penalizedPath.fit[[iter_path]]@unpenalized,
                        penalizedPath.fit[[iter_path]]@penalized,
                        penalizedPath.fit[[iter_path]]@nuisance$sigma2)



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

#### data are not normogonalized
penalized.fit <- penalized(response = df.data[,all.vars(formula.lvm)[1]], 
                           penalized = df.data[,all.vars(formula.lvm)[-1]], 
                           steps = "Park")
unlist(lapply(penalized.fit, function(x){x@lambda1}))

##
p1lvm.fit <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fix.sigma = TRUE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data))
p1lvm.fit$opt$message

p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 621.94685, fix.sigma = TRUE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data))
coef(p1lvm.fit)

##
p1lvm.fit <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fix.sigma = FALSE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data))
p1lvm.fit$opt$message

p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 10, fix.sigma = FALSE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data))
coef(p1lvm.fit)

#### check regularization path - fix.sigma ####

p1lvm.fit <- estimate(plvm.model,  data = df.data_norm, regularizationPath = TRUE, fix.sigma = FALSE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data_norm))
p1lvm.fit$opt$message

p1lvm.fit <- estimate(plvm.model,  data = df.data_norm, lambda1 = p1lvm.fit$opt$message$lambda[11], fix.sigma = FALSE,
                      control = list(constrain = TRUE, iter.max = 1000, data = df.data_norm))
coef(p1lvm.fit)

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



##
penalized.fit <- penalized(response = df.data2[,"Y1"], 
                           penalized = df.data2[,c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")], 
                           steps = "Park")
unlist(lapply(penalized.fit, function(x){x@lambda1}))

penalized.fit <- penalized(response = df.data2[,"Y1"], 
                           penalized = df.data2[,c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")], 
                           lambda1 =  27.77252)


lvm.simple <- lvm(Y1 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)
plvm.simple <- penalize(lvm.simple)

simple.fit <- estimate(plvm.simple,  data = df.data2, lambda1 = 27.77252)
coef(simple.fit)


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
