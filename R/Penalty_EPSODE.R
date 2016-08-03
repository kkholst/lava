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
#' @param resolution_lambda1 the first value is the maximum relative difference in parameter between two steps. 
#' If this lead to a too small step, the second value is used as the minimum change in penalization parameter between two steps.
#' @param nstep_max the maximum number of iterations
#' @param ode.method the type of method to use to solve the ode (see the documentation of deSolve:::ode)
#' @param lambda2 L2 penalization parameter
#' @param group.lambda1 group of lambda to be penalized together (!! not functional now !!)
#' @param control additional options to be passed to the proximal algorithm
#' 
#' @references 
#' Zhou 2014 - A generic Path Algorithm for Regularized Statistical Estimation


EPSODE <- function(beta, beta_lambda0, beta_lambdaMax, objective, gradient, hessian, V, lambda2, group.lambda1, 
                   indexPenalty, indexNuisance, 
                   resolution_lambda1, increasing, stopLambda, stopParam,
                   nstep_max = min(length(beta)*50,1e4), 
                   ode.method = "euler", control, tol.0 = 1e-8, trace){
  
  #### preparation
  n.coef <- length(beta)
  lambda2_save <- lambda2
  lambda2 <- rep(0, n.coef)
  lambda2[indexPenalty] <- lambda2_save 
  envir <- environment()
  
  ## lambda
  res <- initLambda_EPSODE(increasing = increasing,
                           gradient = gradient, beta = beta_lambdaMax, indexPenalty = indexPenalty, indexNuisance = indexNuisance)
  seq_lambda1 <- res$seq_lambda
  stepLambda1 <- res$stepLambda
  if(increasing == FALSE){resolution_lambda1 <- -resolution_lambda1}
  if(length(resolution_lambda1) != 2){
    stop("EPSODE : argument \'resolution_lambda1\' must have length 2 \n",
         "proposed length: ",length(resolution_lambda1),"\n")}
  
  ## reestimate beta
  if(any(lambda2 > 0)){
    
    proxOperator <- function(x, step){  
      control$proxOperator(x, step = step,
                           lambda1 = seq_lambda1, lambda2 = lambda2, test.penalty1 = group.lambda1, test.penalty2 = lambda2>0, 
                           index.constrain = indexNuisance, type.constrain = control$constrain, expX = control$proxGrad$expX)
    }
   
    # may not be ok - check whether constrains are needed 
    beta <- do.call("proxGrad",
                    list(start = beta, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                         step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta,  force.descent = control$proxGrad$force.descent, trace = FALSE, 
                         iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, method = control$proxGrad$method))$par
    
  }
  
  ## constrain 
  if(length(indexNuisance) > 0){
    res <- initSigmaConstrain(beta, constrain = control$constrain, indexNuisance = indexNuisance)
    beta <- res$start
    indexAllCoef <- res$indexAllCoef
  }else{
    constrain <- NULL
    indexAllCoef <- 1:n.coef
  }
  
  #### initialization
  iter <- 1
  test.ncv <- TRUE
  
  ## res
  M.beta <- rbind(beta)
  setNE <- intersect(which(V %*% beta < -tol.0),  indexPenalty)
  setZE <- intersect(which(abs(V %*% beta) < tol.0), indexPenalty)
  setPE <- intersect(which(V %*% beta > tol.0), indexPenalty)
  seq_index <- NA
  
  if(trace>=0){
    cat("Penalisation path using the EPSODE algorithm \n", sep = "")
    if(length(indexNuisance) > 0){
      cat(" * fixed coef : \"",paste(setdiff(names(beta), names(beta)[indexAllCoef]), collapse = "\" \""),"\" \n", sep = "")
      cat(" * value coef : ",paste(beta[setdiff(names(beta), names(beta)[indexAllCoef])], collapse = " ")," \n", sep = "")
    }
    pb <- utils::txtProgressBar(min = 0, max = length(indexPenalty), initial = 0, style = 3)
  }
  
  
  #### main loop
  while(iter < nstep_max && test.ncv>0){
   
    ## current parameters
    iterLambda1 <- seq_lambda1[iter]
    if(iter > 1 && iterLambda1 == seq_lambda1[iter-1]){
      iterLambda1 <- iterLambda1 + stepLambda1*tail(resolution_lambda1,1)
    }
    iterBeta <- M.beta[iter,]
    
    if(length(setZE)>0){
      iterBeta[setZE] <- 0
    }
    #### estimate uz and Uz
    uz <- rep(0, n.coef)
    if(length(setNE)>0){uz <- uz - colSums(V[setNE,,drop = FALSE])}
    if(length(setPE)>0){uz <- uz  + colSums(V[setPE,,drop = FALSE])}
    
    Uz <- V[setZE,indexAllCoef,drop = FALSE]
    Uz_pen <- V[setZE,indexPenalty,drop = FALSE]# Uz[,indexAllCoef %in% indexPenalty, drop = FALSE]
    
    if(length(Uz_pen)>0){ # in case of ill conditionned problem
      B <- pracma::nullspace(Uz)
      iUz_pen <- solve(Uz_pen %*% t(Uz_pen)) %*% Uz_pen # MASS::ginv(Uz_pen) # or (Uz_pen t(Uz_pen))^-1 Uz_pen
    }else{
      B <- NULL
      iUz_pen <- NULL
    }
    
    ## Solve ODE 
    lambda.ode <- seq_len(2000)
    cv.ode <- NULL
    bridge.ode <- rbind(c(iter = 0, step = 0, lambda = iterLambda1, iterBeta))
    
    # H1 <- hessian(iterBeta)
    # H2 <- hessianGaussianO(iterBeta)
    # Hdiff <- H1 - H2
    # attr(Hdiff,"grad") <- attr(H1,"grad") - attr(H2,"grad")
    # print(Hdiff)
    
    res.error <- try(deSolve::ode(y = iterBeta, 
                                  times = lambda.ode, 
                                  func = EPSODE_odeBeta, method = ode.method,
                                  parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, 
                                              lambda2 = lambda2, indexPenalty = indexPenalty, indexAllCoef = indexAllCoef,
                                              uz = uz, Uz = Uz, Uz_pen = Uz_pen, iUz_pen = iUz_pen, B = B, 
                                              resolution = resolution_lambda1, envir = envir)
    ), silent = TRUE)
    # browser()
    
    ## second chance in case of multiple events
    # if(!is.null(cv.ode) && cv.ode$cv["cv.sign"]>1){
    #   bridge.odeS <- bridge.ode
    #   iterBeta2 <- bridge.ode[nrow(bridge.odeS)-2,-(1:3)]
    #   lambda.ode2 <- bridge.ode[nrow(bridge.odeS)-2,3]
    #   
    #   lambda.ode <- seq_len(10000)
    #   cv.ode <- NULL
    #   bridge.ode <- rbind(c(iter = 0, step = 0, lambda = lambda.ode2, iterBeta2))
    #   
    #   res.error <- try(deSolve::ode(y = iterBeta2, 
    #                                 times = lambda.ode, 
    #                                 func = EPSODE_odeBeta, method = ode.method,
    #                                 parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, 
    #                                             lambda2 = lambda2, indexPenalty = indexPenalty, indexAllCoef = indexAllCoef,
    #                                             uz = uz, Uz = Uz, Uz_pen = Uz_pen, iUz_pen = iUz_pen, B = B, 
    #                                             resolution = c(resolution_lambda1[1],resolution_lambda1[2]/10), envir = envir)
    #   ), silent = TRUE)
    # }
    
   #  cat("\n iteration ",iter,"\n")
   # print(cv.ode)
    
    ## update 
    if(is.null(cv.ode)){
      seq_index <- c(seq_index, NA)
    }else if(cv.ode$cv["cv.sign"]==1){
      setNE <- setdiff(setNE,  cv.ode$index)
      setPE <- setdiff(setPE, cv.ode$index)
      setZE <- union(setZE, cv.ode$index)
      seq_index <- c(seq_index, cv.ode$index)
    }else if(cv.ode$cv["cv.constrain"]){
      setZE <- setdiff(setZE, cv.ode$index)
      if(cv.ode$cv["s"]>0){
        setPE <- union(setPE, cv.ode$index) 
      }else{
        setNE <- union(setNE, cv.ode$index)  
      }
      seq_index <- c(seq_index, cv.ode$index)
    }else {
        stop(res.error[1])
    }
    
    newLambda1 <- unname(bridge.ode[nrow(bridge.ode),3])
    M.beta <- rbind(M.beta, bridge.ode[nrow(bridge.ode),-(1:3)])
    seq_lambda1 <- c(seq_lambda1, newLambda1)
    iter <- iter + 1
    
    ## cv
    if(stepLambda1 > 0){
      test.ncv <- (length(setNE) > 0 || length(setPE) > 0 )
      if(!is.null(stopLambda) && newLambda1>=stopLambda){test.ncv <- -1}
      if(!is.null(stopParam) && length(setZE)>=stopParam){test.ncv <- -1}
    }else {
      test.ncv <- (newLambda1!=0) && length(setZE)>0
      if(!is.null(stopLambda) && newLambda1<=stopLambda){test.ncv <- -1}
      if(!is.null(stopParam) && (length(setNE)+length(setPE))>=stopParam){test.ncv <- -1}
    }
    if(trace>=0){utils::setTxtProgressBar(pb, value = if(increasing){length(setZE)}else{length(setNE)+length(setPE)})}
    
  }
  if(trace>=0){close(pb)}
 
  #### post treatment
  if(iter >= nstep_max && trace>=0){
    warning("EPSODE algorithm: maximum number of steps reached \n")
  }
  if(increasing == FALSE && !is.null(beta_lambda0) && test.ncv==0){
    M.beta <- rbind(M.beta,
                    unname(beta_lambda0))
    seq_lambda1 <- c(seq_lambda1, 0)
    seq_index <- c(seq_index, NA)
  }
  
  #### export
  seq_lambda1 <- unname(seq_lambda1)
  rownames(M.beta) <- NULL
  df <- as.data.frame(cbind(lambda1.abs = if(length(indexNuisance) == 0){NA}else{seq_lambda1}, 
                            lambda1 = if(length(indexNuisance) == 0){seq_lambda1}else{NA}, 
                            lambda2.abs = lambda2_save, 
                            lambda2 = NA,
                            indexChange = unname(seq_index),
                            M.beta))
  return(df)
}

EPSODE_odeBeta <- function(t, y, ls.args){
  
  bridge <- get("bridge.ode", envir = ls.args$envir)
  lambda <- bridge[tail(which(bridge[,1]==(t-1)),1),3]
 
  #### test knot
  index <- NULL
  if(length(ls.args$setNE)>0){index <- c(index, ls.args$setNE[which(y[ls.args$setNE] > 0)])} ## any negative parameter that becomes positive: stop algorithm
  if(length(ls.args$setPE)>0){index <- c(index, ls.args$setPE[which(y[ls.args$setPE] < 0)])} ## any positive parameter that becomes negative: stop algorithm
  
  if(length(index) == 1){ 
    assign("cv.ode", 
           value = list(param = y, lambda = lambda, index = index, cv = c(cv.sign = TRUE, cv.constrain = FALSE, s = NA)), 
           envir = ls.args$envir)
    assign("bridge.ode",
           value =  rbind(bridge,c(t, NA, lambda, y)),
           envir = ls.args$envir)
    stop("cv \n")
  }else if(length(index) > 1){
    assign("cv.ode", 
           value = list(param = y, lambda = lambda, index = NA, cv = c(cv.sign = length(index), cv.constrain = FALSE, s = NA)), 
           envir = ls.args$envir)
    assign("bridge.ode",
           value =  rbind(bridge,c(t, NA, NA, y)),
           envir = ls.args$envir)
    stop("EPSODE_odeBeta: multiple events \n",
         "increase resolution \n")
  }
  
  #### Hessian and gradient
  Hfull <- ls.args$hessian(y)
  H <- Hfull[ls.args$indexAllCoef,ls.args$indexAllCoef, drop = FALSE]
  attr(H, "grad")  <- attr(Hfull, "grad")[ls.args$indexAllCoef, drop = FALSE]
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
    H_m1 <- try(solve(H), silent = TRUE)
    if(is.matrix(H_m1)){ # 
      # all coef
      H_m1 <- solve(H)
      R <- solve(ls.args$Uz %*% H_m1 %*% t(ls.args$Uz))
      Q <- H_m1 %*% t(ls.args$Uz) %*% R
      P <- H_m1 - Q %*% ls.args$Uz %*% H_m1 
      
      # only penalized coef
      H_m1 <- solve(H[ls.args$indexAllCoef %in% ls.args$indexPenalty, ls.args$indexAllCoef %in% ls.args$indexPenalty,drop = FALSE])
      R <- solve(ls.args$Uz_pen %*% H_m1 %*% t(ls.args$Uz_pen))
      Q <- H_m1 %*% t(ls.args$Uz_pen) %*% R 
      
    }else{ # singular H matrix
      if(all(H<1e-12)){
        assign("cv.ode", 
               value = list(param = y, lambda = lambda, index = NA, cv = c(cv.sign = FALSE, cv.constrain = FALSE, s = NA)), 
               envir = ls.args$envir)
        assign("bridge.ode",
               value =  rbind(bridge,c(t, NA, NA, y)),
               envir = ls.args$envir)
        stop("EPSODE_odeBeta: all values in the hessian are below 1e-12 \n",
             "try using a larger step \n")
      }
      BHB <-  t(ls.args$B) %*% H %*% ls.args$B
      P <- ls.args$B %*% solve(BHB) %*% t(ls.args$B)
      Q <- ls.args$iUz_pen
      
    }
    G <- attr(H, "grad")[ls.args$indexAllCoef %in% ls.args$indexPenalty, drop = FALSE]
    
    ## check constrains
    s <- - t(Q) %*% ( (1 / lambda) * G + ls.args$uz[ls.args$indexPenalty,drop = FALSE])
    
    if(t > 1 && any( abs(s) > 1)){
      index <- which.max(abs(s))
      assign("cv.ode",
             value =  list(param = y, lambda = lambda, index = ls.args$setZE[index], cv = c(cv.sign = FALSE, cv.constrain = TRUE, s = s[index])),
             envir = ls.args$envir)
      assign("bridge.ode",
             value =  rbind(bridge,c(t, NA, lambda, y)),
             envir = ls.args$envir)
      stop("cv \n")
    }
    
  }
  
  ## export
  Puz <- rep(0, length(y))
  if(length(ls.args$setZE)==0){
    Puz[ls.args$indexAllCoef] <- P %*% ls.args$uz[ls.args$indexAllCoef, drop = FALSE]
  }else{
    Puz[setdiff(ls.args$indexAllCoef,ls.args$setZE)] <- P[-ls.args$setZE,-ls.args$setZE, drop = FALSE] %*% ls.args$uz[setdiff(ls.args$indexAllCoef,ls.args$setZE), drop = FALSE]  
  }
  
  if(ls.args$resolution[1]>0){
    rdiff.max <- max(abs(Puz[setdiff(1:length(y), ls.args$setZE)]/y[setdiff(1:length(y), ls.args$setZE)]))
    normTempo <- max(ls.args$resolution[1]/rdiff.max,ls.args$resolution[2])
  }else{
    lPlus <- (- t(Q) %*% G)/(1 + t(Q) %*% ls.args$uz[ls.args$indexPenalty,drop = FALSE])
    lMinus <- (- t(Q) %*% G)/(-1 + t(Q) %*% ls.args$uz[ls.args$indexPenalty,drop = FALSE])
    lAll <- c(lPlus[lPlus<lambda],lMinus[lMinus<lambda])
    nextKnot <- lAll[which.min(lambda-lAll)]
    normTempo <- min( (lambda-nextKnot)*ls.args$resolution[1], ls.args$resolution[2])
  }
  
  assign("bridge.ode",
         value =  rbind(bridge,c(t, normTempo, lambda+normTempo, y)),
         envir = ls.args$envir)
  # if(any(is.na(Puz)) || any(abs(Puz)>1)){print(Puz);print(y); browser()}
  
  return(list(-Puz*normTempo))
}



initLambda_EPSODE <- function(increasing, gradient, beta, indexPenalty, indexNuisance){
  
  
  if(increasing){
    seq_lambda <- 0 
    
    stepLambda <- max( abs(-gradient(beta) )[indexPenalty] )
    if(length(indexNuisance) > 0){stepLambda <- stepLambda * beta[indexNuisance[1]]}
    
    
  }else{
    
    seq_lambda <-  max( abs(-gradient(beta) )[indexPenalty] ) * 1.1 # initialisation with the fully penalized solution
    if(length(indexNuisance) > 0){seq_lambda <- seq_lambda * beta[indexNuisance[1]]}
    
    stepLambda <- -seq_lambda
    
  }
  
  return(list(seq_lambda = seq_lambda,
              stepLambda = stepLambda))
}



