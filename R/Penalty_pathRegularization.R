LassoPath_lvm <- function(beta0, hessianLv, gradientLv, 
                          indexPenalty, indexNuisance, fix.nuisance,
                          gradientPen = NULL, iter_max = length(beta0)*10,
                          dataY, dataX){
  
  #### orthogonalize data
  res <- normData(dataX = dataX)
  dataX <- res$data
  Xcenter <- res$Xcenter
  
  scalePen <- rep(1, length(beta0))
  names(scalePen) <- names(beta0)
  Xcenter_param <- rep(NA, length(beta0))
  names(Xcenter_param) <- names(beta0)
  
  for(iter_x in 1:length(res$sd.X)){
    index <- grep(pattern = names(res$sd.X)[iter_x], names(scalePen))  
    if(length(index)>0){
      scalePen[index] <- res$sd.X[iter_x]
      Xcenter_param[index] <- Xcenter[iter_x]
    }
  }
  scalePen[indexNuisance] <- 0
  Xcenter_param[indexNuisance] <- 0
  
  
  #### derivative of the likelihood
  n.data <- NROW(dataX)
  #   S <- cov(cbind(dataY,dataX))*(n.data-1)/n.data
  #   mu <- apply(cbind(dataY,dataX),2,mean)
  
  # gradientLv_save <- gradientLv
  
  if(missing(hessianLv)){
    hessianLv <- function(coef, ...){
      
      XX <- cbind(1,dataX)
      YY <- dataY
      
      sigma2 <- coef[length(coef)]
      beta <- coef[-length(coef)]
      residuals <- YY-XX %*% cbind(beta)
      
      hess_sigma2 <- + n/(2*sigma2^2) - 1/(sigma2^3) * crossprod(residuals, residuals)
      hess_sigma2FDbeta <- - 1/(sigma2^2) * crossprod(XX, residuals)
      hess_beta <- - 1/(sigma2) * crossprod(XX, XX)
      
      return(cbind( rbind(hess_beta,t(hess_sigma2FDbeta)),
                    rbind(hess_sigma2FDbeta,hess_sigma2)
      ))
    }
  }
  
  if(missing(gradientLv)){
    gradientLv <- function(coef, ...){
      
      XX <- cbind(1,dataX)
      YY <- dataY
      
      beta <- coef[-length(coef)]
      residuals <- YY-XX %*% cbind(beta)
      sigma2 <- coef[length(coef)]
      
      grad_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * crossprod(residuals, residuals)
      grad_beta <- + 1/(sigma2) * crossprod(XX, residuals)
      
      return(drop(c(grad_beta, grad_sigma2)))
    }
  }
  #   objectiveLv2 <- function(coef, ...){
  #     
  #     XX <- cbind(1,dataX)
  #     YY <- dataY
  #     
  #     beta <- coef[-length(coef)]
  #     residuals <- YY-XX %*% cbind(beta)
  #     sigma2 <- coef[length(coef)]
  #     
  #     sum(dnorm(residuals, mean = 0, sd = sqrt(sigma2), log = TRUE))
  #     
  #     grad_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * crossprod(residuals, residuals)
  #     grad_beta <- + 1/(sigma2) * crossprod(XX, residuals)
  #     
  #     return(drop(c(grad_beta, grad_sigma2)))
  #   }
  
  # coef <- M.beta[1,]
  # gradientLv2 <- function(coef, ...){
  
  #### See Klaus thesis formula 21 page 30 
  #     XX <- cbind(1,dataX)
  #     YY <- dataY
  #     S_XY <- cov(cbind(XX, YY))
  #     mu_XY <- apply(cbind(YY, XX), 2, mean)
  
  #     beta <- coef[-length(coef)]
  #     residuals <- mu_XY[1]- mu_XY[-1] %*% cbind(beta)
  #     sigma2 <- coef[length(coef)]
  #     
  #     grad_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * crossprod(residuals, residuals)
  #     grad_beta <- + 1/(sigma2) * crossprod(XX, residuals)
  
  #     return(drop(c(grad_beta, grad_sigma2)))
  #   }
  
  #   gradientLv_save(beta0)
  #   gradientLv(beta0)
  
  
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
    M.beta[1,indexNuisance] <- 1  
  }
  M.beta[1,1] <- mean(dataY)
  
  V.lambda <- max( abs( gradientLv(M.beta[1,]))[indexPenalty]  * scalePen[indexPenalty] )
  
  cv <- FALSE
  iter <- 1
  
  while(iter < iter_max && cv == FALSE){
    lambda_tempo <- V.lambda[iter]*1/scalePen
    lambda_tempo[-indexPenalty] <- 0
    
    resNode <- nextNode_lvm(hessianLv = hessianLv, gradientLv = gradientLv,  gradientPen = gradientPen,
                            beta = M.beta[iter,], lambda1 = lambda_tempo, 
                            indexPenalty = indexPenalty, indexNuisance = indexNuisance)
    
    newLambda <- V.lambda[iter] * (1 - resNode$gamma)
    
    if(newLambda < 0 || is.infinite(newLambda)){
      cv <- TRUE
    }else{
    V.lambda <- c(V.lambda, newLambda)
    M.beta <- rbind(M.beta, resNode$beta)
    iter <- iter + 1
    }
  }
  
  #### restore the original scale 
  M.beta[,-indexNuisance] <- sweep(M.beta[,-indexNuisance], MARGIN = 2, FUN = "/", STATS = scalePen[-indexNuisance])
  index.col <- which(!is.na(Xcenter_param))
  M.beta[,-index.col] <- M.beta[,-index.col] - Xcenter_param[setdiff(index.col, indexNuisance)] %*% t(M.beta[,setdiff(index.col, indexNuisance)])
  
  #### export
  return(data.frame(lambda = V.lambda, 
                    setNames(as.data.frame(M.beta), c("Intercept", colnames(dataX) ))
  ))
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
  
  grad_B.A <- solve(hess_Lv[set_A,set_A], grad_Pen[set_A]) * lambda1[set_A] 
  grad_Rho <- hess_Lv[,set_A] %*% grad_B.A / lambda1 
  
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


normData <- function(data, formula, dataX = NULL){
  
  if(!missing(data) && !missing(formula)){
    dataX <- as.matrix(data[,all.vars(formula.lvm)[-1]])  
  }
  
  ones <- rep(1,nrow(dataX))
  Xcenter <- solve(crossprod(ones), crossprod(ones, dataX))
  dataX <- dataX - cbind(ones) %*% Xcenter
  sd.X <- sqrt(apply(dataX, 2, var)*(nrow(dataX)-1)/nrow(dataX))
  dataX <- sweep(dataX, MARGIN = 2, FUN = "/", 
                 STATS = sd.X)
  
  
  if(!missing(data) && !missing(formula)){
    dataX <- data.frame(data[,all.vars(formula.lvm)[1], drop = FALSE],dataX)
  }
  
  return(list(data = dataX,
              sd.X = sd.X,
              Xcenter = Xcenter))
}