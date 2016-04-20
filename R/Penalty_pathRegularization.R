LassoPath_lvm <- function(beta0, objectiveLv, hessianLv, gradientLv, 
                          indexPenalty, indexNuisance,
                          sd.X, base.lambda1, lambda2, fix.nuisance, proxOperator, control,
                          gradientPen = NULL, iter_max = length(beta0)*2){

  p <- length(beta0)
  
  if(is.null(gradientPen)){
    gradientPen <- function(beta, grad_Lv, ...){
      ifelse(beta == 0, -sign(grad_Lv), -sign(beta))  # because grad_Lv in lava is -grad(Lv)
    }
  }
  
  #### initialisation
  M.beta <- matrix(0, nrow = 1, ncol = p)
  colnames(M.beta) <- names(beta0)
  M.beta[1,] <- beta0
  if(fix.nuisance == TRUE){
    if(control$constrain){
      M.beta[1,indexNuisance] <- 0
    }else{
      M.beta[1,indexNuisance] <- 1  
    }
  }
  
  V.lambda <- max( abs(-gradientLv(M.beta[1,]) * sd.X)[indexPenalty] )

  cv <- FALSE
  iter <- 1
  
  #### main loop
  while(iter < iter_max && cv == FALSE){
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
      
    if(any(lambda2>0)){ ## correction of the beta due to the L2 correction
      resNode$beta <- do.call("ISTA",
                              list(start = resNode$beta, proxOperator = proxOperator, hessian = hessianLv, gradient = gradientLv, objective = objectiveLv,
                                   lambda1 = newLambda*base.lambda1, lambda2 = lambda2, constrain = resNode$beta[indexNuisance],
                                   iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))$par
    }
      
    if(fix.nuisance == FALSE){ ## estimation of sigma
      warperObj <- function(sigma2){
        pp <- resNode$beta
        pp[indexNuisance] <- sigma2
        return(objectiveLv(pp))
      }
      
      warperGrad <- function(sigma2){
        pp <- resNode$beta
        pp[indexNuisance] <- sigma2
        return(sum(abs(gradientLv(pp)[indexNuisance])))
      }
      
#       test1 <- optim(par = resNode$beta[indexNuisance], 
#                      fn = warperObj,
#                      gr = warperGrad,
#                      method = "Brent", lower = 0, upper = M.beta[1,indexNuisance])$par
      
      suppressWarnings(
      resNode$beta[indexNuisance] <- optim(par = resNode$beta[indexNuisance], 
                                           fn = warperGrad, 
                                           lower = rep(1e-12,length(indexNuisance)), upper = M.beta[1, indexNuisance])$par
      )
    
      ## update lambda given that the product lambda*sigma must be a constant
      newLambda <- median(newLambda * M.beta[iter,names(resNode$beta[indexNuisance])] / resNode$beta[indexNuisance])
    }
    
    ## correction step
    resNode$beta <- do.call("ISTA",
                            list(start = resNode$beta, proxOperator = proxOperator, hessian = hessianLv, gradient = gradientLv, objective = objectiveLv,
                                 lambda1 = newLambda*base.lambda1, lambda2 = lambda2, constrain = if(fix.nuisance){resNode$beta[indexNuisance]}else{NULL},
                                 iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))$par
    
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
  grad_Lv <- -gradientLv(beta)  # because grad_Lv in lava is -grad(Lv)
  hess_Lv <- -hessianLv(beta)  # because hess_Lv in lava is -hess(Lv)
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

EPSODE <- function(gradient, hessian, beta, V, lambda,
                   indexPenalty){
  
  n.beta <- length(beta)
  setNE <- intersect(which(V %*% beta < 0),  indexPenalty)
  setZE <- intersect(which(V %*% beta == 0), indexPenalty)
  setPE <- intersect(which(V %*% beta > 0), indexPenalty)
  
  ODEpath <- function(t, y, ls.args){
    H <- ls.args$hessian(y)
    H_m1 <- solve(H)
    Uz <- ls.args$V[ls.args$setZE,]
    UzHm1Uz_m1 <- solve(Uz %*% H_m1 %*% t(Uz))
    uz <- - colSums(ls.args$V[ls.args$setNE,]) +  colSums(ls.args$V[ls.args$setPE,])
    
    P <- H_m1 - H_m1 %*% t(Uz) %*% UzHm1Uz_m1 %*% Uz %*% H_m1 
    return(list(- P %*% uz))
  }
  
  iter <- 1
  seq_lambda <- lambda
  M.beta <- rbind(beta)
  
  #### iterations
  # while(iter < 1000 && (length(setNE) > 0 || length(setPE) )){}
  iterLambda_m1 <- 0
  iterLambda <- seq_lambda[iter,]
  iterBeta <- M.beta[iter,-1]
  
  ## Solve ODE
  res.ode <- ode(y = iterBeta, times = c(iterLambda_m1,iterLambda), func = ODEpath,  # hmax
                 list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE)
  )
  
  iterBeta <- res.ode[1,-1]
  
  M.beta <- rbind(M.beta,
                  iterBeta)
                  
  ## find next node
  H <- hessian(iterBeta)
  H_m1 <- solve(H)
  G <- gradient(iterBeta)
  Uz <- V[setZE,]
  uz <- - colSums(V[setNE,]) +  colSums(V[setPE,])
  Q <- H_m1 %*% t(Uz) %*% solve(Uz %*% H_m1 %*% t(Uz))
  
  lambda_tempo1 <- ( - Q[indexPenalty,] %*% G[indexPenalty] ) / (1 + Q[indexPenalty,] %*% uz[indexPenalty] )
  lambda_tempo1 <- lambda_tempo1[lambda_tempo1>=0]
  lambda_tempo2 <- ( - Q[indexPenalty,] %*% G[indexPenalty] ) / (- 1 + Q[indexPenalty,] %*% uz[indexPenalty] )
  lambda_tempo2 <- lambda_tempo2[lambda_tempo2>=0]
  
  iterLambda <- min(lambda_tempo1,lambda_tempo2)
  
  if(iterLambda < 0 || is.infinite(iterLambda)){
    cv <- TRUE
  }else{
    seq_lambda <- c(seq_lambda, iterLambda)
  }
  
  ## udpate sets
  setNE <- intersect(which(V %*% beta < 0),  indexPenalty)
  setZE <- intersect(which(V %*% beta == 0), indexPenalty)
  setPE <- intersect(which(V %*% beta > 0), indexPenalty)
  
  ####
  iter <- iter + 1
}
