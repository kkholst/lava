pathRegularization <- function(X, Y, start, objective, gradient, hessian, index.penalized, ...){
  
  index.varCoef <- grep(",", names(start), fixed = TRUE)
  n.coefmean <- ncol(X)
  lambda1 <- rep(0, length(start))
  dots <- list(...)
  
  ## max penalty 
  seq_lambda <- max(abs(gradient(start, penalty = NULL)))
  lambda1[index.penalized] <- seq_lambda[1]
  A <- 1
  I <- setdiff(1:n.coefmean, A)
  
  startPath <- do.call(dots$control$proxGrad.method,
                       list(start = start, step = dots$step, proxOperator = dots$proxOperator, gradient = gradient, 
                            lambda1 = lambda1, lambda2 = dots$lambda2,
                            iter.max = dots$control$iter.max, abs.tol = dots$control$abs.tol, rel.tol = dots$control$rel.tol))$par
  
  ## first non zero parameter
  mu <- X %*% startPath[1:n.coefmean]
  XY <- t(X) %*% (Y - mu)
  max_XY <- max(abs(XY)[-A])
  newA <- which.min(abs(max_XY - abs(XY)))
  A <- c(A, newA)
  I <- setdiff(1:n.coefmean, A)
  
  knot1 <- startPath
  seq_lambda <- c(abs(gradient(knot1, penalty = NULL)[newA]),seq_lambda)
  iter <- 0
  
#   while(knot1[newA] == 0 && iter < 10){
    lambda1[index.penalized] <- seq_lambda[1] - iter * 0.999
   
    knot1 <- do.call(dots$control$proxGrad.method,
                     list(start = startPath, step = dots$step, proxOperator = dots$proxOperator, gradient = gradient, 
                          lambda1 = lambda1, lambda2 = dots$lambda2,
                          iter.max = dots$control$iter.max, abs.tol = dots$control$abs.tol, rel.tol = dots$control$rel.tol))$par
    
    iter <- iter + 1
  # }
  
  mu <- X %*% knot1[1:n.coefmean]
  XY <- t(X) %*% (Y - mu)
  max_XY <- max(abs(XY)[-A])
  newA <- which.min(abs(max_XY - abs(XY)))
  
  hessian2 <- hessian
#   hessian2 <- function(x, penalty){return(diag(1,length(x)))}
#   hessian2 <- function(x, penalty){numDeriv::hessian(func = objective, x = x)}
  newlambda <- calcLambda(coef = knot1, gradient = gradient, hessian = hessian2, index.penalized = index.penalized, A = setdiff(A, 1), newA = newA)
  seq_lambda <- c(newlambda, seq_lambda)
  lambda1[index.penalized] <- seq_lambda[1] 
  
  knot2 <- do.call(dots$control$proxGrad.method,
                   list(start = knot1, step = dots$step, proxOperator = dots$proxOperator, gradient = gradient, 
                        lambda1 = lambda1, lambda2 = dots$lambda2,
                        iter.max = dots$control$iter.max, abs.tol = dots$control$abs.tol, rel.tol = dots$control$rel.tol))$par
  
#   lambda1[penalty$index.coef]  <- abs(gradient(startPath, penalty = NULL)[A]) - 5.3448 # 72.86485 # 78.20965
 
  cat("regularization path: ",paste(seq_lambda, collapse = " ")," \n")
  cat("A set              : ",paste(A, collapse = " ")," \n")
  cat("coef at each knot  : ",paste(knot2, collapse = " ")," \n")
  cat("                   : ",paste(knot1, collapse = " ")," \n")
  cat("                   : ",paste(startPath, collapse = " ")," \n")
  
}


# pathLasso <- function(coef, n.coef, X, A){
#   powMath <- function(x, n){with(eigen(x), vectors %*% (values^n * t(vectors)))}
#   
#   XA <- X[,A,drop = FALSE]
#   OneA <- matrix(1,ncol = 1, nrow = length(A))
#   GA <- t(XA) %*% XA 
#   GAm1 <- solve(GA)
#   AA <- as.numeric(powMath(t(OneA) %*% GAm1 %*% OneA, -1/2))
#   wA <- AA * GAm1 %*% OneA
#   
#   ## LARS
#   a <- t(X) %*% uA
#     
#   index.candidates <- setdiff(1:n.pcoef, A)
#   candidates <- c( (C - cHat[index.candidates])/(AA - a[index.candidates]), 
#                    (C + cHat[index.candidates])/(AA + a[index.candidates]) )
#   gamma_hat <- min(candidates[candidates>0])
#   
#   ## Lasso
#   d <- rep(0,n.coef)
#   d[A] <- wA*sign(coef[A])
#   gamma_tilde <- -coef/d
#   
#   return(list(gamma_hat = gamma_hat,
#               gamma_tilde = gamma_tilde))
# }

calcLambda <- function(coef, gradient, hessian, index.penalized, A, newA){
  
  Hbeta <- as.vector(solve(hessian(coef, penalty = NULL)) %*% sign(coef))
  Hbeta[-A] <- 0
  
  fctTempo <- function(x){
    Gbetah <- abs(gradient(coef + x * Hbeta, penalty = NULL))
    Gbetah_A <- median(Gbetah[A])
    if( all( Gbetah_A > Gbetah[setdiff(index.penalized,A)] )  ){
      diffG <- Gbetah[newA] - Gbetah_A
      return(abs(diffG))
    }else{
      return(10^5)
    }
  }
  
  gamma <- optim(par = 0, fn = fctTempo, lower = 0,
                 method = "L-BFGS-B")$par
  
  return(median(abs(gradient(coef + gamma * Hbeta, penalty = NULL)[A])))
  
}

pathLasso <- function(Y, X, beta, Gbeta, A, I){
  
  ## 2
  gammaA <- -(beta[A]/Gbeta[A])
  
  ## 3
  ymXb <- Y - X %*% beta
  XGbeta <- X %*% Gbeta
  gammaI <- apply(expand.grid(i = A, j = I), 1,
                  function(x){
                    
                    numerator <- c( t(X[,x[1]] + X[,x[2]]) %*% ymXb,
                                    t(X[,x[1]] - X[,x[2]]) %*% ymXb
                    )
                    
                    denominator <- c( t(X[,x[1]] + X[,x[2]]) %*% XGbeta,
                                      t(X[,x[1]] - X[,x[2]]) %*% XGbeta
                    )
                    
                    return(numerator/denominator)
                  }
  )
  
  ## 4
  gamma <- min(c(gammaA,gammaI)[c(gammaA,gammaI)>0])
  
  ## 5-10
  index.addI <- which(abs(gamma-gammaA) < 1e-12)
  index.addA <- which(abs(gamma-gammaI) < 1e-12)
  
  if( length(index.addI) + length(index.addA) > 1){
    stop("pathLasso: several simulatneous update \n")
  }else if( length(index.addA)==1 ){
    A <- c(A, index.addA)
    I <- setdiff(I, index.addA)
  } else if( length(index.addI)==1 ){
    I <- c(I, index.addI)
    A <- setdiff(A, index.addI)
  }else{
    stop("pathLasso: no update \n")
  }
  
  ## 11
  beta <- beta + gamma * Gbeta
  
  ## 12
  Gbeta[A] = - solve(2 * t(X[,A,drop = FALSE]) %*% X[,A,drop = FALSE]) * sign(beta[A])
  
  return(list(beta = beta,
              Gbeta = Gbeta,
              gamma = gamma,
              I = I,
              A = A))
}



# pathLars <- function(Y, X, mu, A){
#   n.pcoef <- ncol(Xall)
#   powMath <- function(x, n){with(eigen(x), vectors %*% (values^n * t(vectors)))}
#   
#   ## initialisation 
#   if(is.null(A)){
#     cHat <- t(X) %*% (Y-mu)
#     C <- max(abs(cHat))
#     A <- which(cHat > C - 1e-12)  
#     
#     ## following steps
#   }else{
#     
#     XA <- X[,A,drop = FALSE]
#     OneA <- matrix(1,ncol = 1, nrow = length(A))
#     GA <- t(XA) %*% XA 
#     GAm1 <- solve(GA)
#     AA <- as.numeric(powMath(t(OneA) %*% GAm1 %*% OneA, -1/2))
#     wA <- AA * GAm1 %*% OneA
#     uA <- XA %*% wA
#     a <- t(X) %*% uA
#     
#     index.candidates <- setdiff(1:n.pcoef, A)
#     candidates <- c( (C - cHat[index.candidates])/(AA - a[index.candidates]), 
#                      (C + cHat[index.candidates])/(AA + a[index.candidates]) )
#     gamma <- min(candidates[candidates>0])
#     A <- c(A,
#            rep(index.candidates,2)[candidates>0][which( abs(candidates - gamma) < 1e-12)]
#     )
#   }
#   
#   return(list(A = A,
#               gamma = gamma))
# }