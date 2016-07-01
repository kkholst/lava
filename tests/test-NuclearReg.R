#### 0- load packages ####
library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(lava)
library(deSolve)


path.lava <- "C:/Users/hpl802/Documents/GitHub/lava" #"C:/Users/hpl802/Downloads/lava-penalization/lava-penalization"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})
source(file.path(path.lava,"tests","FCT.R"))

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

#### 1- Setting ###
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

#res <- simForm(n.obs, xmax, ymax, radius = 0.1, distance = "canberra", coords.centered = FALSE)
# fields:::image.plot(1:xmax, 1:ymax, res$X)
# fields:::image.plot(1:xmax, 1:ymax, res$distCenter+1)

n.confounder <- 5
gamma <- rep(1, n.confounder)
Z <- matrix(rnorm(n.obs*n.confounder), nrow = n.obs, ncol = n.confounder)
Y <- Z %*% gamma + X %*% betaI + rnorm(n.obs)
beta <- setNames(c(0,gamma,betaI,1), c("Y",paste0("Z",1:n.confounder),paste0("X",1:n.coord),"Y,Y"))
formula.lvm <- as.formula( paste0("Y~", paste0("X",1:n.coord, collapse = "+") ) )
# lvm.model <- lvm(formula.lvm)
df.data <- data.frame(Y=Y,data.frame(Z),data.frame(X))

#eplvm.model <- estimate(plvm.model, df.data, lambda1 = lambda1)
objectiveLV(beta)

lambda1 <- 1e3
lambda1.vec <- setNames(rep(0,length(beta)), names(beta))
lambda1.vec[paste0("X",1:n.coord)] <- lambda1
test.penalty1 <- names(lambda1.vec) %in% paste0("X",1:n.coord)
beta[test.penalty1>0] <- 0

proxOperator <-  function(x, step){
  x[test.penalty1] <- proxNuclear(x = x[test.penalty1], step = step,
                                  lambda = unique(lambda1.vec[test.penalty1]), 
                                  nrow = xmax, ncol = ymax)
  return(x)
}


resLV <- proxGrad(start = beta, proxOperator,
                hessian = hessianLV, gradient = gradientLV, objective = objectiveLV,
                step = 1, BT.n = 200, BT.eta = 0.5, force.descent = FALSE, constrain = setNames(c(1), c("Y,Y")),#setNames(c(0,1), c("Y","Y,Y")), 
                iter.max = 50, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

B.LS <- matrix(coef(glmnet:::glmnet(x = X, y = Y, family = "gaussian", alpha = 0, lambda = 0), s = 0)[-1],
               nrow = xmax, ncol = ymax, byrow = TRUE)
B.lambda <- matrix(resLV$par[test.penalty1>0], nrow = xmax, ncol = ymax, byrow = TRUE)


test_that("LVM vs pLVM with lasso", {
  expect_equal(object=dfNuclear(B.LS = B.LS, B.lambda = B.LS, lambda = 0, lambda2 = lambda2),
               expected=rep(0,length(coef(eplvm.fit_tempo1))),
               tolerance=0.001,scale=1)    
})


# image(B.LS)
# image(B.lambda)


lambda2 <- 0
B.LS <- matrix(coef(glmnet:::glmnet(x = X, y = Y, family = "gaussian", alpha = 0, lambda = lambda2), s = 0)[-1],
               nrow = xmax, ncol = ymax, byrow = TRUE)


dfNuclear(B.LS = B.LS, B.lambda = B.lambda, lambda = lambda1, lambda2 = lambda2)





## Not run: 
# fitting SNP-BLUP, i.e. a ridge regression on all the markers across the genome
#
SNP.BLUP.result <- bigRR(y = y, X = X, Z = scale(Z)[,1:100])


resMC <- proxGrad(start = beta, proxOperator,
                hessian = hessianMC, gradient = gradientMC, objective = objectiveMC,
                step = 1, BT.n = 200, BT.eta = 0.5, force.descent = TRUE, constrain = NULL,
                iter.max = 50, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

resLV$par[test.penalty1>0]
resLV$par[1:6]
resMC$par[1:6]
# library(MRIaggr)
# multiplot(coords, resLV$par[test.penalty1>0])
# multiplot(coords, resMC$par[test.penalty1>0])

proxOperatorLasso <-  function(x, step){
  mapply(proxL1, x = x, step = step, lambda = lambda1.vec, test.penalty = test.penalty1)
}



resLasso <- proxGrad(start = beta, proxOperatorLasso,
                hessian = hessianO, gradient = gradientO, objective = objectiveO,
                step = 1, BT.n = 200, BT.eta = 0.8, constrain = setNames(c(0,1), c("Y","Y,Y")), 
                iter.max = 100, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)

# multiplot(coords, resLasso$par[test.penalty1>0])


#### lasso example
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
                        hessian = hessianO, gradient = gradientO, objective = objectiveO,
                        step = 1, BT.n = 20, BT.eta = 0.8, constrain = NULL, 
                        iter.max = 100, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = FALSE)


coef(eplvm.model) - resExternal$par
