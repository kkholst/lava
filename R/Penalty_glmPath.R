#### OVERVIEW
# glmPath: Compute the regularization path for a penalized regression 
# nextNode_lvm: Find the next value of the regularization parameter

#' @title  Compute the regularization path for a penalized regression 
#' 
glmPath <- function(beta0, objectiveLv, hessianLv, gradientLv, 
                    indexPenalty, indexNuisance,
                    sd.X, base.lambda1, lambda2, group.lambda1,
                    control, iter.max = length(beta0)*2){
                          
  p <- length(beta0)
  
  gradientPen <- function(beta, grad_Lv, ...){
    ifelse(beta == 0, -sign(grad_Lv), -sign(beta))  # because grad_Lv in lava is -grad(Lv)
  }
 
  #### initialisation
  M.beta <- matrix(0, nrow = 1, ncol = p)
  colnames(M.beta) <- names(beta0)
  M.beta[1,] <- beta0
  if(!is.null(control$proxGrad$fixSigma)){
    if(control$constrain == TRUE){
      M.beta[1,indexNuisance] <- 0
    }else{
      M.beta[1,indexNuisance] <- 1  
    }
  }
  
  V.lambda <- max( abs(-gradientLv(M.beta[1,]) * sd.X)[indexPenalty] )
  
  cv <- FALSE
  iter <- 1
  
  #### main loop
  while(iter < iter.max && cv == FALSE){
    cat("*")

     ## prediction step: next breakpoint assuming linear path and constant sigma
    resNode <- nextNode_lvm(hessianLv = hessianLv, gradientLv = gradientLv,  gradientPen = gradientPen,
                            beta = M.beta[iter,], lambda1 = V.lambda[iter] * base.lambda1,
                            indexPenalty = indexPenalty, indexNuisance = indexNuisance)
    
    newLambda <- V.lambda[iter] * (1 - resNode$gamma)
    
    if(newLambda < 0 || is.infinite(newLambda)){
      cv <- TRUE
    }else{
      # cat(newLambda," : ",paste(resNode$beta, collapse = " - "),"\n")
      
       if(any(lambda2>0)){ ## update beta value with L2 penalization
        resNode$beta <- do.call("ISTA",
                                list(start = resNode$beta, proxOperator = control$proxOperator, hessian = hessianLv, gradient = gradientLv, objective = objectiveLv,
                                     lambda1 = newLambda*base.lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, constrain = control$proxGrad$fixSigma,
                                     step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, trace = FALSE, 
                                     iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$proxGrad$fast))$par
      }
    
       
      if(is.null(control$proxGrad$fixSigma)){ ## estimation of sigma
        resNode$beta[indexNuisance] <- optim.nuisance(coef = resNode$beta, 
                                                      index.sigma2 = indexNuisance,
                                                      sigma_max = control$proxGrad$sigmaMax,
                                                      gradient = gradientLv)
      
        ## update lambda given that the product lambda*sigma must be a constant
        newLambda <- median(newLambda * M.beta[iter,names(resNode$beta[indexNuisance])] / resNode$beta[indexNuisance])
      }
      
      ## correction step
      if(control$regPath$correctionStep && (all(lambda2 == 0) || is.null(control$proxGrad$fixSigma))){
        
      resNode$beta <- do.call("ISTA",
                              list(start = resNode$beta, proxOperator = control$proxOperator, hessian = hessianLv, gradient = gradientLv, objective = objectiveLv,
                                   lambda1 = newLambda*base.lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, constrain = control$proxGrad$fixSigma,
                                   step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, trace = FALSE, 
                                   iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$proxGrad$fast))$par
      }
      
      ## update
      V.lambda <- c(V.lambda, newLambda)
      M.beta <- rbind(M.beta, resNode$beta)
    } 
    
    iter <- iter + 1
  }
  cat("\n")
  
  #### export
  return(data.frame(lambda1 = V.lambda, lambda2 = NA, M.beta))
}

#' @title  Find the next value of the regularization parameter
nextNode_lvm <- function(hessianLv, gradientLv, gradientPen,
                         beta, lambda1, indexPenalty, indexNuisance){
  
  ##
  hess_Lv <- -hessianLv(beta)  # because hess_Lv in lava is -hess(Lv)
  grad_Lv <- -attr(hess_Lv, "grad")  # because grad_Lv in lava is -grad(Lv) // gradientLv(beta)
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