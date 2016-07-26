#### objective / gradient / hessian
objectiveLV <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  lv <-  - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(epsilon) %*% epsilon
  
  return(-as.numeric(lv))
}

gradientLV <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  Xint <- cbind(1,X)
  sigma2 <- coef[length(coef)]
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * t(epsilon) %*% epsilon    
  gradient_beta <- + 1/(sigma2) * t(Xint) %*% epsilon  
  
  return(-as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianLV <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  G <- gradientLV(coef = coef, Y = Y, X = X)
  
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
  attr(H,"grad") <- G
  
  
  return(-H)
}

objectiveMC <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  mc <-  t(epsilon) %*% epsilon
  
  return(as.numeric(mc))
}

gradientMC <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  gradient_sigma2 <- 0
  gradient_beta <- - t(Xint) %*% epsilon  
  
  return(as.numeric(c(gradient_beta,gradient_sigma2)))
}

hessianLv <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  n <- length(Y)
  G <- gradientMC(coef = coef, Y = Y, X = X)
  
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  hessian_beta <- t(Xint) %*% Xint
  
  hessian_sigma2 <- 0
  hessian_sigma2FDbeta <- rep(0, length(beta))
  H <- cbind( rbind(hessian_beta,t(hessian_sigma2FDbeta)),
              rbind(hessian_sigma2FDbeta,hessian_sigma2)
  )
  attr(H,"grad") <- G
  
  return(H)
}

coef2.penalized <- function(x, iter_lambda){
  
  if(is.list(x)){
    
    if(!missing(iter_lambda)){
      x <- x[[iter_lambda]]
    }else{
      res <- lapply(x, function(model){
        c(model@lambda1,
          model@lambda2,
          model@unpenalized, 
          model@penalized, 
          model@nuisance$sigma2)
      })
      
      Mres <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
      colnames(Mres) <- c("lambda1","lambda2",names(x[[1]]@unpenalized),names(x[[1]]@penalized),"sigma2")
      return(Mres)
    }
    
  } 
  
  coef <- c(x@unpenalized, 
            x@penalized, 
            x@nuisance$sigma2)
  n.coef <- length(coef)
  names.response <- all.vars(x@formula$penalized)[1]
  names(coef)[1] <- names.response
  names(coef)[2:(n.coef-1)] <- paste(names.response, names(coef)[2:(n.coef-1)], sep = "~")
  names(coef)[length(names(coef))] <- paste(names.response, names.response, sep = ",")
  return(coef)
}


validLVM <- function(x, GS = NULL){
  library(penalized)
  
  if(is.null(GS)){
    lambda1 <- x$penalty$lambda1
    lambda1 <- lambda1*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]
    lambda2 <- x$penalty$lambda2
    lambda2 <- lambda2*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]
    
    GS <- penalized(as.formula(paste0(endogenous(x),"~.")), 
                              lambda1 = lambda1, 
                              lambda2 = lambda2, 
                              data = x$data$model.frame, trace = FALSE)
  }
  
  diffCoef <- coef(x) - coef2.penalized(GS)
  
  return(diffCoef)
}