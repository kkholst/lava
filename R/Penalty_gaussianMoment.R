penalized_objective.lvm <- lava:::gaussian_objective.lvm
penalized_gradient.lvm <- lava:::gaussian_gradient.lvm
penalized_hessian.lvm <- lava:::gaussian_hessian.lvm

penalized1_objective.lvm <- lava:::gaussian1_objective.lvm
penalized1_gradient.lvm <- lava:::gaussian1_gradient.lvm
penalized1_hessian.lvm <- lava:::gaussian1_hessian.lvm
# gaussian1_gradient.lvm <- function(x,p,...) {
#   myg <- function(p1) gaussian_objective.lvm(x,p=p1,...)
#   numDeriv::jacobian(myg,p)
# }
# gaussian1_hessian.lvm <- function(x,p,...) {
#   myg <- function(p1) gaussian_objective.lvm(x,p=p1,...)
#   numDeriv::hessian(myg,p)
# }

penalized2_objective.lvm <- lava:::gaussian2_objective.lvm
penalized2_gradient.lvm <- lava:::gaussian2_gradient.lvm
penalized2_hessian.lvm <- lava:::gaussian2_hessian.lvm

penalized3_objective.lvm <- lava:::gaussian3_objective.lvm
penalized3_gradient.lvm <- lava:::gaussian3_gradient.lvm
penalized3_hessian.lvm <- lava:::gaussian3_hessian.lvm

penalized_method.lvm <- "proxGrad"#lava:::gaussian_method.lvm # nlminb2
penalized_logLik.lvm <- lava:::gaussian_logLik.lvm


#### additional

lvGaussian <- function(coef, Y, X, var = NULL){
  
  if(!is.null(var)){coef <- c(coef, var)}

  n <- length(Y)
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
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


