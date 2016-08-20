gaussian_numHessian.lvm <- function(x,p,method,...) {
  
  myg <- function(p1) lava:::gaussian_gradient.lvm(x,p=p1,...)
  I <- numDeriv::jacobian(myg,p, method = method)
  return( (I+t(I))/2 )
  # same as
  # dots <- list(...)
  # dots$type <- "hessian"
  # return(do.call(lava:::gaussian_hessian.lvm, c(list(x,p),dots)))
  
  # system.time(
  #   H0 <- hessianGaussianO(p)
  # )
  # system.time(
  #   H1 <- numDeriv::jacobian(myg,p)
  # )
  # system.time(
  #   H2 <- numDeriv::hessian(myo,p)
  # )
  # print(c(range(H0-H1), " | ",range(H0-H2)))
  # print(c(range(H0-H1)))
  # return(H1)
}

penalized_objective.lvm <- lava:::gaussian_objective.lvm # function(...){gaussian_dlv.lvm(derivative = 0, ...)}#
penalized_gradient.lvm <- lava:::gaussian_gradient.lvm # function(...){gaussian_dlv.lvm(derivative = 1, ...)}#
penalized_hessian.lvm <-  lava:::gaussian_hessian.lvm#gaussian_numHessian.lvm#function(...){gaussian_dlv.lvm(derivative = 2, ...)}#

penalized1_objective.lvm <- lava:::gaussian1_objective.lvm
penalized1_gradient.lvm <- lava:::gaussian1_gradient.lvm
penalized1_hessian.lvm <-  penalized1_hessian.lvm <- function(...){gaussian_dlv.lvm(derivative = 2, ...)}

penalized2_objective.lvm <- lava:::gaussian2_objective.lvm
penalized2_gradient.lvm <- lava:::gaussian2_gradient.lvm
penalized2_hessian.lvm <- lava:::gaussian2_hessian.lvm

penalized3_objective.lvm <- lava:::gaussian3_objective.lvm
penalized3_gradient.lvm <- lava:::gaussian3_gradient.lvm
penalized3_hessian.lvm <- lava:::gaussian3_hessian.lvm



numDeriveSimple_objective.lvm <- lava:::gaussian_objective.lvm
numDeriveSimple_gradient.lvm <- lava:::gaussian_gradient.lvm
numDeriveSimple_hessian.lvm <- function(...){gaussian_numHessian.lvm(method = "simple", ...)}

numDeriveRichardson_objective.lvm <- lava:::gaussian_objective.lvm
numDeriveRichardson_gradient.lvm <- lava:::gaussian_gradient.lvm
numDeriveRichardson_hessian.lvm <- function(...){gaussian_numHessian.lvm(method = "Richardson", ...)}


explicit_objective.lvm <- lava:::gaussian_objective.lvm
explicit_gradient.lvm <- lava:::gaussian_gradient.lvm
explicit_hessian.lvm <-  penalized1_hessian.lvm <- function(...){gaussian_dlv.lvm(derivative = 2, ...)}



penalized_method.lvm <- "proxGrad"
penalized_logLik.lvm <- lava:::gaussian_logLik.lvm

#### additional ####
# here S -> \hat{Sigma}
#      C -> \Omega
gaussian_dlv.lvm <- function(x, p, data, S, mu, n, derivative = 0, ...){
  
  mp <- modelVar(x, p = p, data = data, ...)
  C <- mp$C
  xi <- mp$xi
  iC <- Inverse(C, det = TRUE)
  detC <- attributes(iC)$det
  attr(iC, "det") <- NULL
  u <- mu - xi
  W <- suppressMessages(crossprod(rbind(u)))
  T <- S + W
  
  if(any(1:2 %in% derivative)){
    # do not call the original function as it the computation of the second derivatives are not working
    D <- deriv.lvm2(x, meanpar = attributes(mp)$meanpar, mom = mp,
                    p = p, mu = mu, mean = mean, second = 2 %in% derivative)
    # lava:::deriv.lvm(x, meanpar = attributes(mp)$meanpar, mom = mp, 
    #                       p = p, mu = mu, mean = mean, second = (2 %in% derivative))
    partial_vec.C <- D$dS
    partial_vec.T <- D$dT
    vec.iC <- as.vector(iC)
    vec.iCTiC <- as.vector(iC %*% T %*% iC)
    
  }
  
  ## storage
  res <- vector(mode = "list", length = length(derivative))
  names(res) <- c("lv","score","hessian")[derivative+1]
  
  #### likelihood
  if(0 %in% derivative){
    res$lv <- n/2 * log(detC) + n/2 * tr(T %*% iC)
    # print(res$lv - lava:::gaussian_objective.lvm(x, p, data, S, mu, n, ...))
  }
  
  #### score
  if(1 %in% derivative){
    vec.iC <- as.vector(iC)
    Grad <- crossprod(partial_vec.C, vec.iCTiC - vec.iC) - crossprod(partial_vec.T, vec.iC) 
    res$score <- as.numeric(-n/2 * Grad)
    
    # check derivative
    # calcC <- function(p){as.vector(modelVar(x, p = p, data = data, ...)$C)}
    # range(partial_vec.C -numDeriv::jacobian(func = calcC, x = as.double(p)))
  }
  
  
  #### hessian
  if(2 %in% derivative){
    iCu <- iC %*% t(u)
    iCuuiC <- tcrossprod(iC %*% t(u))
    
    partial_vec.iC <- - (iC %x% iC) %*% partial_vec.C # ok
    partial_vec.xi <- D$dxi # ok
    
    n.rows <- sqrt(nrow(partial_vec.iC))
    ls.tempo <- (lapply(1:n.rows, function(x){partial_vec.iC[x*(1:n.rows),]}))
    partial_iC <- abind:::abind(ls.tempo, along = 3)
    
    partial2_vec.C <- D$d2vecS # D$d2vecS
    partial2_vec.xi <- NULL
    
    B.33 <- (t(iCuuiC) %x% iC) %*%  partial_vec.C
    B.34 <- (iCu %x% iC) %*%  partial_vec.xi
    B.35 <- (iC %x% iCu) %*%  partial_vec.xi
    B.36 <- (iC %x% iCuuiC) %*%  partial_vec.C
    partial_vec.iCTiC <- B.33+B.34+B.35+B.36
    
    Hess.1 <- t(partial_vec.iC) %*% partial_vec.C
    Hess.2 <- t(partial_vec.iCTiC) %*% partial_vec.C
    Hess.3a <- 0 # vec.iC %*% partial2_vec.C 
    Hess.3b <- 0 # vec.iCTiC %*% partial2_vec.C
    Hess.4 <- - 2 * t(partial_vec.xi) %*% iC %*% partial_vec.xi
    Hess.5 <- - 2 * apply(partial_iC, MARGIN = 3, function(x){u %*% x}) %*% partial_vec.xi
    Hess.6 <- 0 #- 2 * t(u) %*% iC %*% partial2_vec.xi
    res$hessian <- -1/2 * (Hess.1 + Hess.2 + Hess.3a + Hess.3b + Hess.4 + Hess.5 + Hess.6)
    browser()
   #### check derivative
   # vec.iC
   test <- FALSE
   if(test){
     calc_iC <- function(p){
       C <- modelVar(x, p = p, data = data, ...)$C
       iC <- Inverse(C, det = TRUE)
       return(as.vector(iC))
     }
     range(partial_vec.iC -numDeriv::jacobian(func = calc_iC, x = as.double(p)))
     
     # vec.xi
     calc_vec.xi <- function(p){
       xi <- modelVar(x, p = p, data = data, ...)$xi
       return(as.vector(xi))
     }
     range(partial_vec.xi -numDeriv::jacobian(func = calc_vec.xi, x = as.double(p)))
   }
  }
  
  if(length(derivative)==1){
    return(res[[1]])
  }else{
    return(res)
  }
}



#### explicit formula ####

lvGaussian <- function(coef, Y, X, var = NULL){
  if(!is.null(var)){coef <- c(coef, var)}
 
  n <- length(Y)
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  if(sigma2<0){return(Inf)}
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  lv <-  - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(epsilon) %*% epsilon
  
  return(-as.numeric(lv))
}

scoreGaussian <- function(coef, Y, X, var = NULL){
  if(!is.null(var)){coef <- c(coef, var)}
  
  n <- length(Y)
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * t(epsilon) %*% epsilon    
  gradient_beta <- + 1/(sigma2) * t(Xint) %*% epsilon  
  
  if(!is.null(var)){gradient_sigma2 <- NULL}
  
  return(-as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianGaussian <- function(coef, Y, X, var = NULL){
  
  G <- scoreGaussian(coef = coef, Y = Y, X = X, var = var)
  
  if(!is.null(var)){coef <- c(coef, var)}
  
  n <- length(Y)
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  hessian_sigma2 <- + n/(2*sigma2^2) - 2/(2*sigma2^3) * t(epsilon) %*% epsilon
  hessian_sigma2FDbeta <- - 1/(sigma2^2) * t(Xint) %*% epsilon
  hessian_beta <- - 1/(sigma2) * t(Xint) %*% Xint
  
  if(!is.null(var)){
    H <- hessian_beta
  }else{
    H <- cbind( rbind(hessian_beta,t(hessian_sigma2FDbeta)),
                c(hessian_sigma2FDbeta,hessian_sigma2)
    )
  }
  attr(H,"grad") <- G
  
  return(-H)
}


