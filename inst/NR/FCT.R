regObj <- function(coef, Y, X){
  
  n <- NROW(Y)
  Xint <- as.matrix(cbind(X))#Xint <- cbind(1,X)
  sigma2 <- as.double(coef[length(coef)])#sigma2 <- coef[1]
  alpha <- NULL#alpha <- coef[2]
  beta <- as.double(coef[-length(coef)])#beta <- coef[-(1:2)]
  epsilon <- Y-Xint %*% cbind(c(alpha,beta))
  obj <-  - n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * t(epsilon) %*% epsilon
  
  return(as.numeric(obj))
}

regFD <- function(coef, Y, X){
  
  n <- NROW(Y)
  Xint <- as.matrix(cbind(X))#Xint <- cbind(1,X)
  sigma2 <- as.double(coef[length(coef)])#sigma2 <- coef[1]
  alpha <- NULL#alpha <- coef[2]
  beta <- as.double(coef[-length(coef)])#beta <- coef[-(1:2)]
  epsilon <- Y-Xint %*% cbind(c(alpha,beta))
  
  FDsigma2 <- - n/(2*sigma2) + 1/(2*sigma2^2) * t(epsilon) %*% epsilon    
  FDbeta <- + 1/(sigma2) * t(Xint) %*% epsilon
  
  return(as.numeric(c(FDbeta,FDsigma2)))
}

regSD <- function(coef, Y, X){
  
  n <- NROW(Y)
  Xint <- as.matrix(cbind(X))#Xint <- cbind(1,X)
  sigma2 <- as.double(coef[length(coef)])#sigma2 <- coef[1]
  alpha <- NULL#alpha <- coef[2]
  beta <- as.double(coef[-length(coef)])#beta <- coef[-(1:2)]
  epsilon <- Y-Xint %*% cbind(c(alpha,beta))
  
  SDsigma2 <- + n/(2*sigma2^2) - 2/(2*sigma2^3) * t(epsilon) %*% epsilon
  FDsigma2FDbeta <- - 1/(sigma2^2) * t(Xint) %*% epsilon
  SDbeta <- - 1/(sigma2) * t(Xint) %*% Xint  
  
  H <- cbind( rbind(SDbeta,t(FDsigma2FDbeta)),
              rbind(FDsigma2FDbeta,SDsigma2)
  )
  
  return(H)
}

