rm(list = ls())

#### 0- load packages ####
library(penalized)
library(data.table)
library(lava)
library(fields)
library(testthat)

library(butils)
package.source("lava", Rcode = TRUE)
source(file.path(butils:::dir.gitHub(),"lava","tests","FCT.R"))



#### 1- Simulation ###
set.seed(10)
n.obs <- 500
xmax <- 64 #10#
ymax <- 64 #10#
center <- c(32,32) #5#
radius <- 10#3#
res <- simForm(n.obs, xmax, ymax, radius = radius)
coords <- res$coords
n.coord <- nrow(coords)
betaI <- as.vector(res$X)
X <- matrix(rnorm(n.obs*n.coord), nrow = n.obs, ncol = n.coord)

# res <- simForm(n.obs, xmax, ymax, radius = 0.1, distance = "canberra", coords.centered = FALSE)
# fields:::image.plot(1:xmax, 1:ymax, res$X)
# fields:::image.plot(1:xmax, 1:ymax, res$distCenter+1)

n.confounder <- 5
gamma <- rep(1, n.confounder)
Z <- matrix(rnorm(n.obs*n.confounder), nrow = n.obs, ncol = n.confounder)
Y <- Z %*% gamma + X %*% betaI + rnorm(n.obs)
formula.lvm <- as.formula( paste0("Y~", paste0("X",1:n.coord, collapse = "+") ) )
# lvm.model <- lvm(formula.lvm)
df.data <- data.frame(Y=Y,data.frame(Z),data.frame(X))
range(unlist(lapply(df.data, sd)))
df.data[] <- scale(df.data)

## display
MRIaggr:::multiplot(as.data.table(coords), betaI) # coeffcient
MRIaggr:::multiplot(as.data.table(coords), X[1,]) # realisation


#### 2- Model ####
names.param <- c("intercept",paste0("Z",1:n.confounder),paste0("X",1:n.coord),"Y,Y")
beta <- setNames(c(0,gamma,betaI,1), names.param)
beta[paste0("X",1:n.coord)] <-  0

test.penaltyMC <- setNames(names(beta) %in% paste0("X",1:n.coord), names(beta))
test.penaltyLV <- setNames(test.penaltyMC[-length(test.penaltyMC)], names(beta)[-length(test.penaltyMC)])


#### glmnet
glmnet.fit <- glmnet:::glmnet(x = cbind(Z,X), y = Y, family = "gaussian", alpha = 0,
                              penalty.factor = test.penaltyLV)

B.LS <- matrix(coef(glmnet.fit, s = 0)[-(1:(n.confounder+1))],
               nrow = xmax, ncol = ymax, byrow = TRUE)
fields:::image.plot(B.LS)

B.LS <- matrix(coef(glmnet.fit, s = 100)[-(1:(n.confounder+1))],
               nrow = xmax, ncol = ymax, byrow = TRUE)
fields:::image.plot(B.LS)


loss.glmnet <- sapply(1:length(glmnet.fit$lambda), function(x){
  sum((glmnet.fit$beta[-(1:n.confounder),x]-betaI)^2)
})
plot(glmnet.fit$lambda,loss.glmnet)


#### estimation using mean square loss 
lambda1 <- 7e2
lambda1.vec <- setNames(rep(lambda1,n.confounder+n.coord+2), names.param)

#res <- objectiveMC(beta)
#res <- gradientMC(beta)
#res <- hessianMC(beta)

proxOperator <-  function(x, step){
  x[test.penaltyMC] <- proxNuclear(x = x[test.penaltyMC], step = step,
                                   lambda = unique(lambda1.vec[test.penaltyMC]), 
                                   nrow = xmax, ncol = ymax)
  return(x)
}


resMC <- proxGrad(start = beta, proxOperator,
                  hessian = hessianMC, gradient = gradientMC, objective = objectiveMC,
                  step = 1, BT.n = 100, BT.eta = 0.8, force.descent = TRUE,
                  iter.max = 50, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.MC <- matrix(resMC$par[test.penaltyMC>0], nrow = xmax, ncol = ymax, byrow = TRUE)
image.plot(B.MC)

#### estimation using lv but at fixed sigma
# res1 <- objectiveMC(beta) - objectiveLV(beta[-length(beta)], var = 1/2)
# res2 <- objectiveMC(beta+1) - objectiveLV(beta[-length(beta)]+1, var = 1/2)
# print(res1-res2)
# res1 <- gradientMC(beta)[-length(beta)] - gradientLV(beta[-length(beta)], var = 1/2)
# res2 <- gradientMC(beta+1)[-length(beta)] - gradientLV(beta[-length(beta)]+1, var = 1/2)
# print(res1-res2)
# res <- hessianMC(beta) - objectiveLV(beta[-length(beta)], var = 1/2)
lambda1 <- 1e2
lambda1.vec <- setNames(rep(lambda1,n.confounder+n.coord+2), names.param)

objectiveNuclear <- function(coef){lvGaussian(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"]), var = 1)}
gradientNuclear <- function(coef){scoreGaussian(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"]), var = 1)}
hessianNuclear <- function(coef){hessianGaussian(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"]), var = 1)}

proxOperator <-  function(x, step){
  x[test.penaltyLV] <- proxNuclear(x = x[test.penaltyLV], step = step,
                                   lambda = unique(lambda1.vec[test.penaltyLV]), 
                                   nrow = xmax, ncol = ymax)
  return(x)
}

resLV <- proxGrad(start = beta[-length(beta)], proxOperator,
                  objective = objectiveNuclear, gradient = gradientNuclear, hessian = hessianNuclear,
                  step = 1, BT.n = 200, BT.eta = 0.5, force.descent = TRUE,
                  iter.max = 10, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.LV <- matrix(resLV$par[test.penaltyLV>0], nrow = xmax, ncol = ymax, byrow = TRUE)
image.plot(B.LV)
image.plot(B.LV-B.MC)


#### estimation using lv
lambda1 <- 1e2
lambda1.vec <- setNames(rep(lambda1,n.confounder+n.coord+2), names.param)

objectiveNuclear <- function(coef){lvGaussian(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"]))}
gradientNuclear <- function(coef){scoreGaussian(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"]))}
hessianNuclear <- function(coef){hessianGaussian(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"]))}

proxOperator <-  function(x, step){
  x[test.penaltyMC] <- proxNuclear(x = x[test.penaltyMC], step = step,
                                   lambda = unique(lambda1.vec[test.penaltyMC]), 
                                   nrow = xmax, ncol = ymax)
  return(x)
}

resLV <- proxGrad(start = beta, proxOperator,
                  objective = objectiveNuclear, gradient = gradientNuclear, hessian = hessianNuclear,
                  step = 1, BT.n = 100, BT.eta = 0.8, force.descent = TRUE,
                  iter.max = 100, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.LV <- matrix(resLV$par[test.penaltyLV>0], nrow = xmax, ncol = ymax, byrow = TRUE)
image.plot(B.LV)
image.plot(B.LV-B.MC)

estimate(lvm.model, lambdaN = 1e3, data = df.data)

#### estimate the degree of freedom to find the best model
dfNuclear(B.LS = B.LS, B.lambda = B.lambda, lambda = lambda1, lambda2 = lambda2)
## then integrate it into lava throught an additional option specifying:
# the contribution of the image term to the log-likelihood, gradient and hessian
# 



# multiplot(coords, resLasso$par[test.penalty1>0])


#### lvm 
test <- FALSE
if(test){
formula.conf <- paste0("Y~",paste0("Z",1:n.confounder,collapse = "+"),"+",paste0("X",1:n.coord,collapse = "+"))
lvm.model <- lvm(as.formula(formula.conf))
penalizeNuclear(lvm.model, coords = coords) <- as.formula(paste0("Y~",paste0("X",1:n.coord,collapse = "+")))

estimate(lvm.model, lambdaN = 1e3, data = df.data)

res <- estimate(lvm.model, lambdaN = 1e3, data = df.data)


}