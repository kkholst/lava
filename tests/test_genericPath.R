rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(deSolve)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)


objectiveO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  lv <-  - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(epsilon) %*% epsilon
  
  return(as.numeric(lv))
}

gradientO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * t(epsilon) %*% epsilon    
  gradient_beta <- + 1/(sigma2) * t(Xint) %*% epsilon  
  
  return(as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianO <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  hessian_sigma2 <- + n/(2*sigma2^2) - 2/(2*sigma2^3) * t(epsilon) %*% epsilon
  hessian_sigma2FDbeta <- - 1/(sigma2^2) * t(Xint) %*% epsilon
  hessian_beta <- - 1/(sigma2) * t(Xint) %*% Xint
  
  H <- cbind( rbind(hessian_beta,t(hessian_sigma2FDbeta)),
              rbind(hessian_sigma2FDbeta,hessian_sigma2)
  )
  
  return(H)
}

####
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seqPark_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seqParkNorm_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))

plvm.RP <- estimate(plvm.model, data = df.data, regularizationPath = TRUE)
plvm.RP$opt$message


### check Q condition


plvm.RP <- estimate(plvm.model, data = df.data, regularizationPath = 1,
                    control = list(start = coef(estimate(lvm.model, data = df.data))))
# plvm.RP$opt$message
plvm.RP <- estimate(plvm.model, data = df.data, regularizationPath = 2,
                    control = list(constrain = FALSE, step_lambda1 = 10, start = coef(estimate(lvm.model, data = df.data))))

plvm.RP <- estimate(plvm.model, data = df.data, regularizationPath = 2,
                    control = list(constrain = FALSE, step_lambda1 = -10, iter.max = 1000))


# 1 86.153347       0 -0.121409118  0.000000e+00  0.000000e+00 0.000000e+00 0.000000e+00 0.0000000000 18.778273
# 2 86.326405       0 -0.121386816  0.000000e+00  0.000000e+00 0.000000e+00 0.000000e+00 0.0003032559 18.738755
# 3 78.974093       0 -0.036645741  0.000000e+00  0.000000e+00 0.000000e+00 4.628701e-09 1.1525888808 12.699449
# 4 77.484400       0  0.007523891  0.000000e+00  0.000000e+00 5.377920e-09 7.754719e-01 1.8912276217  7.826057
# 5  7.398226       0  0.052847809  0.000000e+00  0.000000e+00 9.050240e-01 1.861319e+00 2.9244174862  3.980388
# 6  7.397491       0  0.052848212  0.000000e+00 -6.790435e-06 9.050291e-01 1.861325e+00 2.9244221768  3.980386
# 7  6.286164       0  0.053473845 -1.692839e-08 -1.026122e-02 9.127400e-01 1.870240e+00 2.9317331449  3.976663
# 8  0.000000       0  0.056772636 -5.865776e-02 -7.109498e-02 9.577815e-01 1.921950e+00 2.9777820863  3.963549
###PB!!!!!!!!!!!!!!
plvm.punctual1 <- estimate(plvm.model, data = df.data, lambda1 = 78.974093)
coef(plvm.punctual1)

plvm.punctual2 <- estimate(plvm.model, data = df.data, lambda1 = 78.974093,
                          control = list(start = coef(estimate(lvm.model, data = df.data))))
coef(plvm.punctual2)

plvm.punctual3 <- estimate(plvm.model, data = df.data, lambda1 = 78.974093,
                           control = list(start = newBeta, step = NULL))
coef(plvm.punctual3)

####
beta <- coef(elvm.model)
V <- diag(0, length(beta))
indexPenalty <- 2:6
diag(V)[indexPenalty] <- 1
# res <- EPSODE(gradient, hessian, beta = start, V = V, indexPenalty = 2:6)

n.beta <- length(beta)
  

  
iter <- 1
seq_lambda <- 0
step_lambda <- 10
M.beta <- rbind(beta)
setNE <- intersect(which(V %*% beta < 0),  indexPenalty)
setZE <- intersect(which(V %*% beta == 0), indexPenalty)
setPE <- intersect(which(V %*% beta > 0), indexPenalty)

#### iterations [1] 86.325999 86.326406 78.974095 77.484401 77.484285  7.398227  7.397492  6.286166  0.000000
while(iter < 20 && (length(setNE) > 0 || length(setPE) > 0 )){
  iterLambda <- seq_lambda[iter]
  iterBeta <- M.beta[iter,]
  if(length(setZE)>0){
    iterBeta[setZE] <- 0
  }
  ## Solve ODE # library(deSolve)
  res.ode <- ode(y = iterBeta, times = seq(iterLambda, iterLambda + step_lambda, length.out = 1000), 
                  func = ODEpath,  
                  parm = list(hessian = hessianO, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE)
  )
  # matplot(res.ode)
  
  #### case 1: shrinkage
  # re run if simultaneous changes or too large coefficient
  index_changeSign <- apply(res.ode[,1+c(setNE,setPE)], 2, function(x){which(diff(sign(x))!=0)[1]})
  if(!all(is.na(index_changeSign))){
    index_newLambda <- min(index_changeSign, na.rm = TRUE)
    index_new0 <- c(setNE,setPE)[which.min(index_changeSign)]
    seq_lambda <- c(seq_lambda, res.ode[index_newLambda,1])
    M.beta <- rbind(M.beta, res.ode[index_newLambda,-1])
    
    ## Update active constrains
    setNE <- setdiff(setNE, index_new0)
    setZE <- union(setZE, index_new0)
    setPE <- setdiff(setPE, index_new0)
    
    
  }else{ 
    seq_lambda <- c(seq_lambda, res.ode[nrow(res.ode),1])
    M.beta <- rbind(M.beta, res.ode[nrow(res.ode),-1])
    
  }
  
  iter <- iter + 1
  #### case 2 
 
  ## Update inactive constrains
#   H <- -hessian(iterBeta)
#   H_m1 <- solve(H)
#   G <- -gradient(iterBeta)
#   Uz <- V[setZE,]
#   uz <- - colSums(V[setNE,]) +  colSums(V[setPE,])
#   Q <- H_m1 %*% t(Uz) %*% solve(Uz %*% H_m1 %*% t(Uz))
#   
#   lambda_tempo1 <- ( - Q[indexPenalty,] %*% G[indexPenalty] ) / (1 + Q[indexPenalty,] %*% uz[indexPenalty] )
#   lambda_tempo1 <- lambda_tempo1[lambda_tempo1>=0]
#   lambda_tempo2 <- ( - Q[indexPenalty,] %*% G[indexPenalty] ) / (- 1 + Q[indexPenalty,] %*% uz[indexPenalty] )
#   lambda_tempo2 <- lambda_tempo2[lambda_tempo2>=0]
#   
#   iterLambda <- min(lambda_tempo1,lambda_tempo2)
#   
#   if(iterLambda < 0 || is.infinite(iterLambda)){
#     cv <- TRUE
#   }else{
#     seq_lambda <- c(seq_lambda, iterLambda)
#   }
#   
#   ## udpate sets
#   setNE <- intersect(which(V %*% beta < 0),  indexPenalty)
#   setZE <- intersect(which(V %*% beta == 0), indexPenalty)
#   setPE <- intersect(which(V %*% beta > 0), indexPenalty)
  
  ####
  
}