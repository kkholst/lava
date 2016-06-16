#### OVERVIEW
# EPSODE: perform the generic Path Algorithm. Because of solve(H) not suited to high dimensional data
# EPSODE_odeBeta: compute the righ side of the ODE linking lambda1 and beta

#' @title Perform the generic Path Algorithm for a LVM
#' @param beta initial values for the parameters
#' @param objective likelihood given by lava. Used to adjust the step parameter when using backtracking
#' @param gradient first derivative of the likelihood given by lava. 
#' @param hessian second derivative of the likelihood given by lava. Only used to estimate the step parameter of the algorithm when step = NULL
#' @param V matrix that left multiply beta to define the penalization (identity corresponds to a standard lasso penalty)
#' @param indexPenalty position of the penalised coefficients in beta
#' @param stepLambda1 the range of lambda values investigated at each step
#' @param resolution_lambda1 the number of lambda values that form the grid with range stepLambda1
#' @param nstep_max the maximum number of iterations
#' @param ode.method the type of method to use to solve the ode (see the documentation of deSolve:::ode)
#' @param lambda2 L2 penalization parameter
#' @param group.lambda1 group of lambda to be penalized together (!! not functional now !!)
#' @param control additional options to be passed to the proximal algorithm
#' 
#' @references 
#' Zhou 2014 - A generic Path Algorithm for Regularized Statistical Estimation


EPSODE <- function(beta, beta_lambdaMax, objective, gradient, hessian, V, lambda2, group.lambda1, 
                   indexPenalty, indexNuisance, 
                   stepLambda1, stepIncreasing, resolution_lambda1 = 1000, nstep_max = min(length(beta)*5,Inf), 
                   ode.method = "euler", control, trace){
  
  #### preparation
  n.coef <- length(beta)
  lambda2_save <- lambda2
  lambda2 <- rep(0, n.coef)
  lambda2[indexPenalty] <- lambda2_save 
  envir <- environment()
   
  #### constrain 
  if(length(indexNuisance) > 0){
    indexAllCoef <- setdiff(1:n.coef, indexNuisance[1])
    beta[indexNuisance] <- beta[indexNuisance]/beta[indexNuisance[1]]
    if(control$constrain){
      beta[indexNuisance] <- log(beta[indexNuisance])
    }
  }else{
    indexAllCoef <- 1:n.coef
  }

  #### initialization
  iter <- 1
  if(is.null(stepLambda1)){
    stepLambda1 <- max( abs(-gradient(beta_lambdaMax) )[indexPenalty] )*if(stepIncreasing){1}else{-1}
    if(length(indexNuisance) > 0){stepLambda1 <- stepLambda1 * sum(beta_lambdaMax[indexNuisance[1]])}
  }
  
  if(stepLambda1 > 0){
    seq_lambda1 <- 0  
  }else{
    seq_lambda1 <-  max( abs(-gradient(beta_lambdaMax) )[indexPenalty] ) * sum(beta_lambdaMax[indexNuisance]) * 1.1 # initialisation with the fully penalized solution
  }
 
  if(any(lambda2 > 0) | stepLambda1 < 0){ ### TOFIX for sum sigma
    if(length(indexNuisance) > 0){
      constrain <- setNames( rep(1 - control$constrain, length(indexNuisance) ), names(beta)[indexNuisance]) 
    }else{
      constrain <- NULL
    }
     beta <- do.call("proxGrad",
                    list(start = beta, proxOperator = control$proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                         lambda1 = seq_lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, constrain = constrain,
                         step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, trace = FALSE, 
                         iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, method = control$proxGrad$method))$par
  }
  
  ##
  M.beta <- rbind(beta)
  setNE <- intersect(which(V %*% beta < 0),  indexPenalty)
  setZE <- intersect(which(V %*% beta == 0), indexPenalty)
  setPE <- intersect(which(V %*% beta > 0), indexPenalty)
  test.ncv <- TRUE
  
  if(trace){
    cat("fixed coef      : \"",paste(setdiff(names(beta), names(beta)[indexAllCoef]), collapse = "\" \""),"\" \n", sep = "")
    cat("value fixed coef: ",paste(beta[setdiff(names(beta), names(beta)[indexAllCoef])], collapse = " ")," \n", sep = "")
  }
   
  #### main loop
  while(iter < nstep_max && test.ncv){
    if(trace){cat("*")}
    
   ## current parameters
    iterLambda1 <- seq_lambda1[iter]
    if(iter > 1 && iterLambda1 == seq_lambda1[iter-1]){
      iterLambda1 <- iterLambda1 + stepLambda1/resolution_lambda1
    }
    iterBeta <- M.beta[iter,]
    
#     cat(iter," ",iterLambda1," ", paste(iterBeta, collapse = " "),"\n")  
    if(length(setZE)>0){
      iterBeta[setZE] <- 0
    }
    
    ## Solve ODE 
    lambda.ode <- seq(iterLambda1, max(0, iterLambda1 + stepLambda1), length.out = resolution_lambda1)
    cv.ODE <- c(cv = FALSE, lambda = lambda.ode[resolution_lambda1], cv.sign = FALSE, cv.constrain = FALSE, s = NA)
   
     res.ode <- deSolve::ode(y = iterBeta, 
                   times = lambda.ode, 
                   func = EPSODE_odeBeta, method = ode.method,
                   parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, 
                               lambda2 = lambda2, indexPenalty = indexPenalty, indexAllCoef = indexAllCoef, 
                               envir = envir)
    )
    index.breakpoint <-  which.min(abs(res.ode[,1] - cv.ODE["lambda"]))
    indexM.breakpoint <- max(1, index.breakpoint - 1)
    indexP.breakpoint <- min(resolution_lambda1, index.breakpoint + 1)
    
    ## more precise estimate
    if(cv.ODE["cv"] == 1){
      lambda.ode <- seq(res.ode[indexM.breakpoint, 1], res.ode[indexP.breakpoint, 1], length.out = resolution_lambda1)
      cv.ODE <- c(cv = FALSE, lambda = lambda.ode[resolution_lambda1], cv.sign = FALSE, cv.constrain = FALSE, s = NA)
      
      res.ode <- deSolve::ode(y = res.ode[indexM.breakpoint,-1], 
                     times = lambda.ode, 
                     func = EPSODE_odeBeta, method = ode.method,
                     parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, 
                                 lambda2 = lambda2, indexPenalty = indexPenalty, indexAllCoef = indexAllCoef, 
                                 envir = envir)
      )
      index.breakpoint <-  which.min(abs(res.ode[,1] - cv.ODE["lambda"]))
      
      if(cv.ODE["cv.sign"] == 1){
        setNE <- setdiff(setNE,  cv.ODE["index"])
        setPE <- setdiff(setPE, cv.ODE["index"])
        setZE <- union(setZE, cv.ODE["index"])
      }else if(cv.ODE["cv.constrain"] == 1){
        setZE <- setdiff(setZE, cv.ODE["index"])
        if(cv.ODE["s"]>0){
          setPE <- union(setPE, cv.ODE["index"]) 
        }else{
          setNE <- union(setNE, cv.ODE["index"])  
        }
      }
      
    }
    #     cat(res.breakpoint$breakpoint.sign, ", ", res.breakpoint$breakpoint.constrain, ": ", res.breakpoint$index_coef,"\n")
    
    ## prepare update
    newLambda1 <- cv.ODE["lambda"]
    newBeta <- res.ode[index.breakpoint,-1]
 
    ## updates
    M.beta <- rbind(M.beta, newBeta)
    seq_lambda1 <- c(seq_lambda1, newLambda1)
    iter <- iter + 1
    if(stepLambda1 > 0){
      test.ncv <- (length(setNE) > 0 || length(setPE) > 0 )
    }else{
      test.ncv <- (newLambda1!=0)#length(setZE)
    }
    
  }
  if(trace){cat("\n")}
  
  #### export
  seq_lambda1 <- unname(seq_lambda1)
  rownames(M.beta) <- NULL
  
  return(as.data.frame(cbind(lambda1.abs = if(length(indexNuisance) == 0){NA}else{seq_lambda1}, 
                             lambda1 = if(length(indexNuisance) == 0){seq_lambda1}else{NA}, 
                             lambda2.abs = lambda2_save, 
                             lambda2 = NA,
                             M.beta)))
}

EPSODE_odeBeta <- function(t, y, ls.args){
  
  #### check cv
  if(get("cv.ODE", envir = ls.args$envir)[1]>0){
    return(list(rep(0,length(y))))
  }
  
  #### test knot
  if(length(ls.args$setNE)>0 && any(y[ls.args$setNE]>0)){ ## any negative parameter that becomes positive: stop algorithm
    index <- ls.args$setNE[which.max(y[ls.args$setNE])]
    assign("cv.ODE", 
           value = c(cv = TRUE, lambda = t, cv.sign = TRUE, cv.constrain = FALSE, index = index, s = NA), 
           envir = ls.args$envir)
    return(list(rep(0,length(y))))
  }
  if(length(ls.args$setPE)>0 && any(y[ls.args$setPE]<0)){ ## any positive parameter that becomes negative: stop algorithm
    index <- ls.args$setPE[which.min(y[ls.args$setPE])]
    assign("cv.ODE", 
           value = c(cv = TRUE, lambda = t, cv.sign = TRUE, cv.constrain = FALSE, index = index, s = NA), 
           envir = ls.args$envir)
    return(list(rep(0,length(y))))
  }
  
  #### estimate uz
  uz <- rep(0, length(y))
  if(length(ls.args$setNE)>0){
    uz <- uz - colSums(ls.args$V[ls.args$setNE,,drop = FALSE])
  }
  if(length(ls.args$setPE)>0){
    uz <- uz  + colSums(ls.args$V[ls.args$setPE,,drop = FALSE])
  }
  
  #### Hessian and gradient
  Hfull <- ls.args$hessian(y)
  H <- Hfull[ls.args$indexAllCoef,ls.args$indexAllCoef, drop = FALSE]
  attr(H, "grad")  <- attr(Hfull, "grad")[ls.args$indexAllCoef, drop = FALSE] ## can we instead set the derivative to 0 ??
  if(any(ls.args$lambda2>0)){
    attr(H, "grad") <- attr(H, "grad") + ls.args$lambda2[ls.args$indexAllCoef, drop = FALSE] * y[ls.args$indexAllCoef, drop = FALSE]
    H[] <- H[] + diag(ls.args$lambda2[ls.args$indexAllCoef, drop = FALSE])
  }
 
  #### estimate Q and P
  
  if(length(ls.args$setZE) == 0){
    R <- NULL
    Q <- NULL
    P <- solve(H)
  }else{
    Uz <- ls.args$V[ls.args$setZE,ls.args$indexAllCoef,drop = FALSE]

### modif high dimension
#     H <- H[-ls.args$setZE,-ls.args$setZE, drop = FALSE]
#     Uz <- Uz[,-ls.args$setZE, drop = FALSE]
    
    H_m1 <- solve(H)
    R <- solve(Uz %*% H_m1 %*% t(Uz))
    Q <- H_m1 %*% t(Uz) %*% R
    P <- H_m1 - Q %*% Uz %*% H_m1 
  }
  
  ## check constrains
  if(!is.null(Q)){
    H_m1 <- solve(H[ls.args$indexAllCoef %in% ls.args$indexPenalty, ls.args$indexAllCoef %in% ls.args$indexPenalty,drop = FALSE])#solve(H[ls.args$indexPenalty, ls.args$indexPenalty,drop = FALSE])  #
    Uz <- Uz[,ls.args$indexAllCoef %in% ls.args$indexPenalty, drop = FALSE]#Uz[,ls.args$indexPenalty, drop = FALSE]                             #
    G <- attr(H, "grad")[ls.args$indexAllCoef %in% ls.args$indexPenalty, drop = FALSE] #attr(H, "grad")[ls.args$indexPenalty, drop = FALSE]                  #
    
    R <- solve(Uz %*% H_m1 %*% t(Uz))
    Q <- H_m1 %*% t(Uz) %*% R  #     Q <- Q[ls.args$indexPenalty,,drop = FALSE]
    s <- - t(Q) %*% ( (1 / t) * G + uz[ls.args$indexPenalty,drop = FALSE])
    
    if(any( abs(s) > 1)){
      index <- which.max(abs(s))
      assign("cv.ODE", 
             value = c(cv = TRUE,  lambda = t, cv.sign = FALSE, cv.constrain = TRUE, index = ls.args$setZE[index], s = s[index]), 
             envir = ls.args$envir)
      return(list(rep(0,length(y))))
    }
    
  }
  
   ## export
  Puz <- rep(0, length(y))
  
  if(length(ls.args$setZE)==0){
    Puz[ls.args$indexAllCoef] <- P %*% uz[ls.args$indexAllCoef, drop = FALSE]
  }else{
    # Puz[setdiff(ls.args$indexAllCoef,ls.args$setZE)] <- P %*% uz[setdiff(ls.args$indexAllCoef,ls.args$setZE), drop = FALSE]  
    Puz[setdiff(ls.args$indexAllCoef,ls.args$setZE)] <- P[-ls.args$setZE,-ls.args$setZE, drop = FALSE] %*% uz[setdiff(ls.args$indexAllCoef,ls.args$setZE), drop = FALSE]  
  }
  
  return(list(-Puz))
}

# EPSODE_odeP <- function(t, y, ls.args){
#   dH <- genD(ls.args$hessian, as.numeric(ls.args$beta))
#   return( (y %x% y) %*% dH %*% y %*% ls.args$uz )
# }

