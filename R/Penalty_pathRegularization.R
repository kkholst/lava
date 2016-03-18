LassoPath_lvm <- function(beta0, hessianLv, gradientLv, indexPenalty, indexNuisance,
                      gradientPen = NULL, iter_max = length(beta0)*10,  scalePen = NULL){
  

#   hessianLv <- function(coef, ...){
#     
#     XX <- cbind(1,dataX_orth)
#     YY <- df.data[,all.vars(formula.lvm)[1]]
#     
#     sigma2 <- coef[length(coef)]
#     beta <- coef[-length(coef)]
#     residuals <- YY-XX %*% cbind(beta)
#     
#     hess_sigma2 <- + n/(2*sigma2^2) - 1/(sigma2^3) * crossprod(residuals, residuals)
#     hess_sigma2FDbeta <- - 1/(sigma2^2) * crossprod(XX, residuals)
#     hess_beta <- - 1/(sigma2) * crossprod(XX, XX)
#     
#     return(cbind( rbind(hess_beta,t(hess_sigma2FDbeta)),
#                   rbind(hess_sigma2FDbeta,hess_sigma2)
#     ))
#   }
  
#   gradientLv <- function(coef, ...){
#     
#     XX <- cbind(1,dataX_orth)
#     YY <- df.data[,all.vars(formula.lvm)[1]]
#     
#     sigma2 <- coef[length(coef)]
#     beta <- coef[-length(coef)]
#     residuals <- YY-XX %*% cbind(beta)
#     
#     grad_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * crossprod(residuals, residuals)
#     grad_beta <- + 1/(sigma2) * crossprod(XX, residuals)
#     
#     return(drop(c(grad_beta, grad_sigma2)))
#   }
  
  p <- length(beta0)
  if(is.null(scalePen)){ scalePen <- rep(1,p)}
    
  if(is.null(gradientPen)){
    gradientPen <- function(beta, grad_Lv, ...){
      ifelse(beta == 0, sign(grad_Lv), sign(beta))  
    }
  }
  
  #### initialisation
  M.beta <- matrix(0, nrow = 1, ncol = p)
  M.gradB <- matrix(0, nrow = 1, ncol = p)
  M.beta[1,] <- beta0
  M.beta[1,indexNuisance] <- 1
  
  V.lambda <- max( abs( gradientLv(M.beta[1,]) * scalePen)[indexPenalty] )
  
  cv <- FALSE
  iter <- 1
  
  while(iter < iter_max && cv == FALSE){
    resNode <- nextNode_lvm(hessianLv = hessianLv, gradientLv = gradientLv,  gradientPen = gradientPen,
                            beta = M.beta[iter,], lambda1 = V.lambda[iter], 
                            scalePen = scalePen, indexPenalty = indexPenalty, indexNuisance = indexNuisance)
    
    newLambda <- V.lambda[iter] * (1 - resNode$gamma)
    
    V.lambda <- c(V.lambda, newLambda)
    M.beta <- rbind(M.beta, resNode$beta)
    
    iter <- iter + 1
    if(newLambda < 0 || is.infinite(newLambda)){
      cv <- TRUE
    }
  }
  
  browser()
   return(data.frame(lambda = V.lambda, setNames(as.data.frame(M.beta), c("Intercept", colnames(dataX) ))))
}

#' @title  Find the next value of the regularization parameter

nextNode_lvm <- function(hessianLv, gradientLv, gradientPen,
                     beta, lambda1, 
                     scalePen, indexPenalty, indexNuisance){

  lambda1 <- lambda1*1/scalePen
  lambda1[-indexPenalty] <- 0
  
  ##
  grad_Lv <- gradientLv(beta)
  hess_Lv <- hessianLv(beta) * scalePen
  grad_Pen <- gradientPen(beta, grad_Lv = grad_Lv)
  
  set_A <- union(setdiff(which(beta!=0),indexNuisance),
                 setdiff(which(abs(grad_Lv)/lambda1 > 1-1e-04),indexNuisance)
  )
  
  grad_B.A <- - solve(hess_Lv[set_A,set_A], grad_Pen[set_A])
  grad_Rho <- hess_Lv[,set_A] %*% grad_B.A
  
  ## find next knot
  gamma1 <- -beta[set_A]/(grad_B.A*lambda1[set_A]*c(0,scalePen)[set_A]) 
  gamma2Moins <- (1 - grad_Lv/lambda1)/(1 + grad_Rho)
  gamma2Moins[c(set_A,indexNuisance)] <- Inf
  gamma2Plus <- (1 + grad_Lv/lambda1)/(1 - grad_Rho)
  gamma2Plus[c(set_A,indexNuisance)] <- Inf
  
  test.min <- any(c(gamma1,gamma2Moins,gamma2Plus)>=0)
  if(test.min){
    gamma <- max(1e-04, min(c(gamma1,gamma2Moins,gamma2Plus)[c(gamma1,gamma2Moins,gamma2Plus)>=0]))
  }else{
    gamma <- Inf
  }
  
  ## update coef
  beta[set_A] <- beta[set_A] + gamma * grad_B.A * lambda1[set_A] * scalePen[set_A]
  
  return(list(gamma = gamma, beta = beta))
}


# pathRegularization <- function(X, Y, start, objective, gradient, hessian, index.penalized, ...){
#   
#   index.varCoef <- grep(",", names(start), fixed = TRUE)
#   n.coefmean <- ncol(X)
#   lambda1 <- rep(0, length(start))
#   dots <- list(...)
#   
#   ## max penalty 
#   seq_lambda <- max(abs(gradient(start, penalty = NULL)))
#   lambda1[index.penalized] <- seq_lambda[1]
#   A <- 1
#   I <- setdiff(1:n.coefmean, A)
#   
#   startPath <- do.call(dots$control$proxGrad.method,
#                        list(start = start, step = dots$step, proxOperator = dots$proxOperator, gradient = gradient, 
#                             lambda1 = lambda1, lambda2 = dots$lambda2,
#                             iter.max = dots$control$iter.max, abs.tol = dots$control$abs.tol, rel.tol = dots$control$rel.tol))$par
#   
#   ## first non zero parameter
#   mu <- X %*% startPath[1:n.coefmean]
#   XY <- t(X) %*% (Y - mu)
#   max_XY <- max(abs(XY)[-A])
#   newA <- which.min(abs(max_XY - abs(XY)))
#   A <- c(A, newA)
#   I <- setdiff(1:n.coefmean, A)
#   
#   knot1 <- startPath
#   seq_lambda <- c(abs(gradient(knot1, penalty = NULL)[newA]),seq_lambda)
#   iter <- 0
#   
#   #   while(knot1[newA] == 0 && iter < 10){
#   lambda1[index.penalized] <- seq_lambda[1] - iter * 0.999
#   
#   knot1 <- do.call(dots$control$proxGrad.method,
#                    list(start = startPath, step = dots$step, proxOperator = dots$proxOperator, gradient = gradient, 
#                         lambda1 = lambda1, lambda2 = dots$lambda2,
#                         iter.max = dots$control$iter.max, abs.tol = dots$control$abs.tol, rel.tol = dots$control$rel.tol))$par
#   
#   iter <- iter + 1
#   # }
#   
#   mu <- X %*% knot1[1:n.coefmean]
#   XY <- t(X) %*% (Y - mu)
#   max_XY <- max(abs(XY)[-A])
#   newA <- which.min(abs(max_XY - abs(XY)))
#   
#   hessian2 <- hessian
#   #   hessian2 <- function(x, penalty){return(diag(1,length(x)))}
#   #   hessian2 <- function(x, penalty){numDeriv::hessian(func = objective, x = x)}
#   newlambda <- calcLambda(coef = knot1, gradient = gradient, hessian = hessian2, index.penalized = index.penalized, A = setdiff(A, 1), newA = newA)
#   seq_lambda <- c(newlambda, seq_lambda)
#   lambda1[index.penalized] <- seq_lambda[1] 
#   
#   knot2 <- do.call(dots$control$proxGrad.method,
#                    list(start = knot1, step = dots$step, proxOperator = dots$proxOperator, gradient = gradient, 
#                         lambda1 = lambda1, lambda2 = dots$lambda2,
#                         iter.max = dots$control$iter.max, abs.tol = dots$control$abs.tol, rel.tol = dots$control$rel.tol))$par
#   
#   #   lambda1[penalty$index.coef]  <- abs(gradient(startPath, penalty = NULL)[A]) - 5.3448 # 72.86485 # 78.20965
#   
#   cat("regularization path: ",paste(seq_lambda, collapse = " ")," \n")
#   cat("A set              : ",paste(A, collapse = " ")," \n")
#   cat("coef at each knot  : ",paste(knot2, collapse = " ")," \n")
#   cat("                   : ",paste(knot1, collapse = " ")," \n")
#   cat("                   : ",paste(startPath, collapse = " ")," \n")
#   
# }