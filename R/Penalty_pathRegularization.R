LassoPath_lvm <- function(beta0, objectiveLv, hessianLv, gradientLv, 
                          indexPenalty, indexNuisance,
                          sd.X, base.lambda, fix.nuisance, exp.nuisance,
                          gradientPen = NULL, iter_max = length(beta0)*10){

  p <- length(beta0)
  
  if(is.null(gradientPen)){
    gradientPen <- function(beta, grad_Lv, ...){
      ifelse(beta == 0, -sign(grad_Lv), -sign(beta))  
    }
  }
  
  #### initialisation
  M.beta <- matrix(0, nrow = 1, ncol = p)
  M.gradB <- matrix(0, nrow = 1, ncol = p)
  M.beta[1,] <- beta0
  if(fix.nuisance == TRUE){
    if(exp.nuisance){
      M.beta[1,indexNuisance] <- 0
    }else{
      M.beta[1,indexNuisance] <- 1  
    }
  }
    
  V.lambda <- max( abs( gradientLv(M.beta[1,]) * sd.X)[indexPenalty] )
  cv <- FALSE
  iter <- 1
  
  while(iter < iter_max && cv == FALSE){
    lambda_tempo <- V.lambda[iter]*base.lambda
    lambda_tempo[-indexPenalty] <- 0
    
    resNode <- nextNode_lvm(hessianLv = hessianLv, gradientLv = gradientLv,  gradientPen = gradientPen,
                            beta = M.beta[iter,], lambda1 = lambda_tempo, 
                            indexPenalty = indexPenalty, indexNuisance = indexNuisance)
    
    newLambda <- V.lambda[iter] * (1 - resNode$gamma)
   
    if(newLambda < 0 || is.infinite(newLambda)){
      cv <- TRUE
    }else{
    V.lambda <- c(V.lambda, newLambda)
    
    if(fix.nuisance == FALSE){
     
#     resNode$beta[indexNuisance] <- nextNuisance_lvm(resNode = resNode, 
#                                                     index = indexNuisance, 
#                                                     objective = objectiveLv, 
#                                                     gradient = gradientLv)
    }
    
    M.beta <- rbind(M.beta, resNode$beta)
    iter <- iter + 1
    }
  }
  
   #### export
  return(data.frame(lambda = V.lambda, M.beta))
}

#' @title  Find the next value of the regularization parameter

nextNode_lvm <- function(hessianLv, gradientLv, gradientPen,
                         beta, lambda1, indexPenalty, indexNuisance){
  
  ##
  grad_Lv <- gradientLv(beta)
  hess_Lv <- hessianLv(beta)
  grad_Pen <- gradientPen(beta, grad_Lv = grad_Lv)
  
  set_A <- union(setdiff(which(beta!=0),indexNuisance),
                 setdiff(which(abs(grad_Lv)/lambda1 > 1-1e-04),indexNuisance)
  )
  
  grad_B.A <- solve(hess_Lv[set_A,set_A, drop = FALSE], grad_Pen[set_A, drop = FALSE] * lambda1[set_A, drop = FALSE])
  grad_Rho <- hess_Lv[,set_A, drop = FALSE] %*% grad_B.A / lambda1 
  
  ## find next knot
  gamma1 <- -beta[set_A]/(grad_B.A*lambda1[set_A]) 
  gamma2Moins <- (1 - grad_Lv/lambda1)/(1 + grad_Rho)
  gamma2Moins[c(set_A,indexNuisance)] <- Inf
  gamma2Plus <- (1 + grad_Lv/lambda1)/(1 - grad_Rho)
  gamma2Plus[c(set_A,indexNuisance)] <- Inf
  
  # cat(paste(gamma1, collpase = " ")," | ",paste(gamma2Moins, collpase = " ")," | ",paste(gamma2Plus, collpase = " ")," | \n \n")
  
  test.min <- any(c(gamma1,gamma2Moins,gamma2Plus)>=0)
  if(test.min){
    gamma <- max(1e-04, min(c(gamma1,gamma2Moins,gamma2Plus)[c(gamma1,gamma2Moins,gamma2Plus)>=0]))
  }else{
    gamma <- Inf
  }
  
  ## update coef
  beta[set_A] <- beta[set_A] + gamma * grad_B.A
  
  return(list(gamma = gamma, beta = beta))
}

nextNuisance_lvm <- function(beta, index, objective, gradient){
 
   fn_warper <- function(x){
    start <- beta
    start[index] <- x
    objective(start)
  }
  gn_warper <- function(x){
    start <- beta
    start[index] <- x
    return(gradient(start)[index])
  }
  
  suppressWarnings(
  resOptim <- optim(par = beta[index],
                    fn = fn_warper, gr = gn_warper, control = list(fnscale = -1))
  )
  
  return(resOptim$par)
}
