#### objective / gradient / hessian
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

hessianMC <- function(coef, Y = df.data$Y, X = as.matrix(df.data[, names(df.data) != "Y"])){
  
  G <- gradientMC(coef = coef, Y = Y, X = X)
  
  n <- length(Y)
  Xint <- cbind(1,X)
  beta <- coef[-length(coef)]
  epsilon <- Y - Xint %*% cbind(beta)
  
  hessian_beta <- t(Xint) %*% Xint
  
  hessian_sigma2 <- 0
  hessian_sigma2FDbeta <- rep(0, length(beta))
  
  H <- cbind( rbind(hessian_beta,t(hessian_sigma2FDbeta)),
              c(hessian_sigma2FDbeta,hessian_sigma2)
  )
  attr(H,"grad") <- G
  
  return(-H)
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
