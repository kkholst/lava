simForm <- function(n.obs, xmax, ymax, radius, center = NULL, coords.centered = TRUE,
         distance = "euclidean"){
  
  if(is.null(center)){
    if(coords.centered == TRUE){
      center <- c(0,0)
    }else{
      center <- c(xmax/2,ymax/2)
    }
  }
  
  coords <- scale(expand.grid(1:xmax, 1:ymax), center = coords.centered, scale = FALSE)
  n.coord <- nrow(coords) 
  
  distCenter <- apply(coords, 1, function(x){dist( rbind(x,center), method = distance)})
  beta <- distCenter<radius
  
  return(list(coords = coords,
              center = center,
              distCenter = matrix(distCenter, nrow = xmax, ncol = ymax),
              X = matrix(beta, nrow = xmax, ncol = ymax)
  ))
}

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

validPath.lvm <- function(x, data, validMean = TRUE, validCov = TRUE, control = list(), ...){
  
  Mall <- getPath(x, getLambda = c("lambda1","lambda2"), rm.duplicated = TRUE)
  Mcoef <- Mall[,-(1:2),drop = FALSE]
  seqLambda1 <- Mall[,"lambda1"]
  seqLambda2 <- Mall[,"lambda2"]
  seqLambda2[is.na(seqLambda2)] <- 0
  n.Lambda <- length(seqLambda1)
  
  if("trace" %in% names(control) == FALSE){
    control$trace <- FALSE
  }
  
  McoefGS <- NULL
  for(iterL in 1:n.Lambda){
    McoefGS <- rbind(McoefGS, coef(estimate(x$x, data = data, control = control,
                                        lambda1 = seqLambda1[iterL], lambda2 = seqLambda2[iterL],
                                        ...)))
  }
  
  return(list(proxAlgo = McoefGS,
              EPSODE = Mcoef,
              diff = Mcoef-McoefGS,
              lambda1 = seqLambda1,
              lambda2 = seqLambda2,
              diff.range = range(Mcoef-McoefGS)))
}

