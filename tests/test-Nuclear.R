#### 0- load packages ####
library(penalized)
library(data.table)
library(lava)
library(fields)
library(testthat)

library(butils)
package.source("lava", Rcode = TRUE)
source(file.path(butils::dir.gitHub(),"lava","tests","FCT.R"))

simForm <- function(n.obs, xmax, ymax, radius, center = NULL, coords.centered = TRUE,
                    distance = "euclidean"){

  if(is.null(center)){
    if(coords.centered == TRUE){
      center <- c(0,0)
    }else{
      center <- c(xmax/2,ymax/2)
    }
  }
  
  coords <- scale(expand.grid(1:xmax, 1:ymax), center = coords.centered, scale = FALSE)
  n.coord <- nrow(coords) 

  distCenter <- apply(coords, 1, function(x){dist( rbind(x,center), method = distance)})
  beta <- distCenter<radius
  
  return(list(coords = coords,
              center = center,
              distCenter = matrix(distCenter, nrow = xmax, ncol = ymax),
              X = matrix(beta, nrow = xmax, ncol = ymax)
              ))
}

#### 1- Simulation ###
n.obs <- 500
xmax <- 64 
ymax <- 64 
center <- c(32,32) 
radius <- 10
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

## display
MRIaggr:::multiplot(as.data.table(coords), betaI) # coeffcient
MRIaggr:::multiplot(as.data.table(coords), X[1,]) # realisation

#### 2- Model ####
names.param <- c("intercept",paste0("Z",1:n.confounder),paste0("X",1:n.coord),"Y,Y")
lambda1 <- 1e3
lambda1.vec <- setNames(rep(lambda1,n.confounder+n.coord+2), names.param)
test.penaltyMC <- names(lambda1.vec) %in% paste0("X",1:n.coord)
test.penaltyLV <- test.penaltyMC[-length(test.penaltyMC)]
beta <- setNames(c(0,gamma,betaI,1), names.param)


#### glmnet
lambda1 <- 1e3
glmnet.fit <- glmnet:::glmnet(x = cbind(Z,X), y = Y, family = "gaussian", alpha = 0,
                              penalty.factor = test.penalty1)

# coef(glmnet.fit)[2:(n.confounder+1),]
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

#beta[test.penalty1>0] <- 0


#### estimation using mean square loss 
res <- objectiveMC(beta)
res <- gradientMC(beta)
res <- hessianMC(beta)

proxOperator <-  function(x, step){
  x[test.penaltyMC] <- proxNuclear(x = x[test.penaltyMC], step = step,
                                   lambda = unique(lambda1.vec[test.penalty1]), 
                                   nrow = xmax, ncol = ymax)
  return(x)
}

resMC <- proxGrad(start = beta, proxOperator,
                  hessian = hessianMC, gradient = gradientMC, objective = objectiveMC,
                  step = 1, BT.n = 200, BT.eta = 0.5, force.descent = FALSE,
                  iter.max = 50, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.MC <- matrix(resMC$par[test.penaltyMC>0], nrow = xmax, ncol = ymax, byrow = TRUE)
image.plot(B.MC)

#### estimation using lv but at fixed sigma
res1 <- objectiveMC(beta) - objectiveLV(beta[-length(beta)], var = 1/2)
res2 <- objectiveMC(beta+1) - objectiveLV(beta[-length(beta)]+1, var = 1/2)
print(res1-res2)
res1 <- gradientMC(beta)[-length(beta)] - gradientLV(beta[-length(beta)], var = 1/2)
res2 <- gradientMC(beta+1)[-length(beta)] - gradientLV(beta[-length(beta)]+1, var = 1/2)
print(res1-res2)
res <- hessianMC(beta) - objectiveLV(beta[-length(beta)], var = 1/2)


objectiveNuclear <- function(coef){objectiveLV(coef,  var = 1)}
gradientNuclear <- function(coef){gradientLV(coef, var = 1)}
hessianNuclear <- function(coef){hessianLV(coef,  var = 1)}

proxOperator <-  function(x, step){
  x[test.penaltyLV] <- proxNuclear(x = x[test.penaltyLV], step = step,
                                   lambda = unique(lambda1.vec[test.penalty1]), 
                                   nrow = xmax, ncol = ymax)
  return(x)
}

resLV <- proxGrad(start = beta[-length(beta)], proxOperator,
                  hessian = hessianNuclear, gradient = gradientNuclear, objective = objectiveNuclear,
                  step = 1, BT.n = 200, BT.eta = 0.5, force.descent = TRUE,
                  iter.max = 10, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.LV <- matrix(resLV$par[test.penaltyLV>0], nrow = xmax, ncol = ymax, byrow = TRUE)
image.plot(B.LV)
image.plot(B.LV-B.MC)



proxOperator <-  function(x, step){
  x[test.penaltyMC] <- proxNuclear(x = x[test.penaltyMC], step = step,
                                   lambda = unique(lambda1.vec[test.penalty1]), 
                                   nrow = xmax, ncol = ymax)
  return(x)
}

resLV <- proxGrad(start = beta, proxOperator,
                  hessian = objectiveLV, gradient = gradientLV, objective = hessianLV,
                  step = 1, BT.n = 200, BT.eta = 0.5, force.descent = TRUE,
                  iter.max = 10, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.LV <- matrix(resLV$par[test.penaltyLV>0], nrow = xmax, ncol = ymax, byrow = TRUE)
image.plot(B.LV)
image.plot(B.LV-B.MC)

#### estimate the degree of freedom to find the best model
dfNuclear(B.LS = B.LS, B.lambda = B.lambda, lambda = lambda1, lambda2 = lambda2)
## then integrate it into lava throught an additional option specifying:
# the contribution of the image term to the log-likelihood, gradient and hessian
# 



# multiplot(coords, resLasso$par[test.penalty1>0])


#### DEBUG - lasso example ####
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <-   as.list( c(rep(0,2),1:3) ) # as.list( c(rep(0,2),0.25,0.5,0.75) ) # 
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)
df.data <- as.data.frame(scale(df.data))

lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

proxOperator <- function(x, step){
  mapply(proxL1, x = x, step = step, lambda = lambda1.vec, test.penalty = test.penalty1)
}
n.coef <- length(coef(plvm.model))

lambda1 <- 10
lambda1.vec <- setNames(rep(0,n.coef), coef(plvm.model))
lambda1.vec[plvm.model$penalty$names.penaltyCoef] <- lambda1
test.penalty1 <- lambda1.vec>0

####
eplvm.model <- estimate(plvm.model, df.data, lambda1 = lambda1, 
                        control = list(start = coef(elvm.model)))
# eplvm.model <- estimate(plvm.model, df.data, lambda1 = lambda1)

resExternal <- proxGrad(start = coef(elvm.model), proxOperator,
                        hessian = hessianLV, gradient = gradientLV, objective = objectiveLV,
                        step = 1, BT.n = 20, BT.eta = 0.8, force.descent = FALSE,
                        iter.max = 100, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = FALSE)


expect_equal(coef(eplvm.model),resExternal$par, tolerance = 1e-7)
