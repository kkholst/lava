proxGrad <- function(start, objective, gradient, hessian,...){
   
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol","proxGrad.method","fix.sigma","test.path")
  dots <- list(...)
  control <- dots$control
  control <- control[names(control) %in% PGcontrols]
  penalty <- dots$penalty
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  n.coef <- length(start)
  
  index.penaltyCoef <- which(names(start) %in% penalty$names.penaltyCoef)
  
  #### definition of the operator
  if(penalty$lambda1 == 0 && penalty$lambda2 == 0){
    proxOperator <- function(x, step, lambda1, lambda2){x}
  }else if(penalty$lambda2 == 0){
    test.penalty <- 1:n.coef %in% index.penaltyCoef
    
    proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty)
    }
  }else if(penalty$lambda1 == 0){
    test.penalty <- 1:n.coef %in% index.penaltyCoef
   
     proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL2, x = x, step = step, lambda = lambda2, test.penalty = test.penalty)
    }
  }else{
    test.penalty1 <- 1:n.coef %in% index.penaltyCoef
    test.penalty2 <- 1:n.coef %in% index.penaltyCoef
    
    proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL2, 
             x = mapply(proxL1, x = x, step = step,  lambda = lambda1, test.penalty = test.penalty1),
             step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
  }
  
     
  #### main
  if(dots$fix.sigma){
    constrain <- rep(1-dots$control$constrain, length(penalty$names.varCoef))
    names(constrain) <- penalty$names.varCoef
    
    if(any(penalty$names.penaltyCoef %in% names(constrain))){
      stop("proxGrad: wrong specification of the penalty \n",
           "cannot penalize the variance parameter when fix.sigma = TRUE \n",
           "requested penalisation on variance parameters: ",penalty$names.penaltyCoef[penalty$names.penaltyCoef %in% names(constrain)],"\n")
    }
    
  }else{
    constrain <- NULL
  }
 
   if( dots$regularizationPath ){
    
    ## path regularization
    resLassoPath <- LassoPath_lvm(beta0 = start, hessianLv = hessian, gradientLv = gradient, objectiveLv = objective,
                                  indexPenalty = index.penaltyCoef, 
                                  indexNuisance = which(names(start) %in% penalty$names.varCoef),
                                  dataY = as.matrix(dots$control$data[,1,drop = FALSE]), dataX = as.matrix(dots$control$data[,-1,drop = FALSE]))
    
    # print(resLassoPath)
    
    res <- list(par = start,
                convergence = 0,
                iterations = 0,
                evaluations = c("function" = 0, "gradient" = 0),
                message = resLassoPath,
                objective = NA)
  }else{
    
    ## one step lasso
    step <- 1/max(abs(eigen(hessian(start, penalty = NULL))$value)) # may be computationnaly expensive 
    
    lambda1 <- rep(0, n.coef)
    lambda1[index.penaltyCoef] <- penalty$lambda1 
    lambda2 <- rep(0, n.coef)
    lambda2[index.penaltyCoef] <- penalty$lambda2 
   
    res <- do.call(dots$proxGrad.method,
                   list(start = start, step = step, proxOperator = proxOperator, gradient = gradient, 
                        lambda1 = lambda1, lambda2 = lambda2, constrain = constrain,
                        iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol))
      
    res$objective <- objective(res$par, penalty = penalty)
  
    if(!is.null(constrain)){
    res2 <- do.call(dots$proxGrad.method,
                   list(start = res$par, step = step, proxOperator = proxOperator, gradient = gradient, 
                        lambda1 = lambda1, lambda2 = lambda2, constrain = res$par[names(res$par) %in% names(constrain) == FALSE],
                        iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol))
    res$par <- res2$par
    }
  }

  ### export
  return(res)
}



ISTA <- function(start, step = NULL, decreaseStep = 0.99, proxOperator, gradient, lambda1, lambda2, constrain,
                 iter.max, abs.tol, rel.tol){
  
  
  ## initialisation
  test.cv <- FALSE
  iter <- 1
  x_k <- start 
  if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){
    x_km1 <- x_k
    step <- step * decreaseStep
    
    x_k <- proxOperator(x = x_km1 - step * gradient(x_km1, penalty = NULL), 
                        step = step, lambda1 = lambda1, lambda2 = lambda2)
    if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
    
    iter <- iter + 1
    absDiff <-  (abs(x_k-x_km1) < abs.tol)
    relDiff <- (abs(x_k-x_km1)/abs(x_k) < rel.tol)
    test.cv <- all(  absDiff + ifelse(is.na(relDiff),0,relDiff) > 0 )
    
  }
  
  # print(x_k - x_km1)
  
  ## export
  message <- if(test.cv){"Sucessful convergence \n"
  }else{
    paste("max absolute/relative difference: ",max(absDiff),"/",max(relDiff)," for parameter ",which.max(absDiff),"/",which.max(relDiff),"\n")
  }
  return(list(par = x_k,
              convergence = as.numeric(test.cv==FALSE),
              iterations = iter,
              evaluations = c("function" = 0, "gradient" = iter),
              message = message
  ))
}

FISTA <- function(start, step = NULL, decreaseStep = 0.99, proxOperator, gradient, lambda1, lambda2, constrain,
                  iter.max, abs.tol, rel.tol){
  
  ## initialisation
  test.cv <- FALSE
  iter <- 1
  x_k <- start 
  t_k <- 1
  if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  
  y_k <- x_k
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){
    x_km1 <- x_k
    t_km1 <- t_k
    step <- step * decreaseStep
  # cat(iter," ")

    x_k <- proxOperator(x = y_k - step * gradient(y_k, penalty = NULL), 
                        step = step, lambda1 = lambda1, lambda2 = lambda2)
    if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  
    t_k <- (1 + sqrt(1 + 4 * t_k^2)) / 2
    y_k <- x_k + (t_km1-1)/t_k * (x_k - x_km1)
    
    iter <- iter + 1
    absDiff <-  abs(x_k-x_km1)
    relDiff <- abs(x_k-x_km1)/abs(x_k)
    test.cv <- all( ( absDiff<abs.tol) + ifelse(is.na(relDiff),0, relDiff < rel.tol) > 0 )
  }
  
  ## export
  message <- if(test.cv){"Sucessful convergence \n"
  }else{
    paste("max absolute/relative difference: ",max(absDiff),"/",max(relDiff, na.rm = TRUE)," for parameter ",which.max(absDiff),"/",which.max(relDiff),"\n")
  }
  return(list(par = x_k,
              convergence = as.numeric(test.cv==FALSE),
              iterations = iter,
              evaluations = c("function" = 0, "gradient" = iter),
              message = message
  ))
}

optimx1 <- function(start, objective, gradient, hessian,...) {
  require(optimx)
  
  optimxcontrols <- c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min","optimx.method")
  dots <- list(...)
  control <- dots$control
  control <- control[names(control)%in%optimxcontrols]
  dots$control <- NULL
  dots$tol.grad_pen <- NULL
  dots$debug <- NULL
  dots$lower <- NULL
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")

  if(!is.null(control$optimx.method)){
    optimx.method <- control$optimx.method
    control$optimx.method <- NULL
   
    optimx.res <- do.call(optimx,
                          c(list(par = start, fn = objective, gr = gradient, method = optimx.method, control = control),
                            dots)
    )
    
  }else{
    control$all.methods <- TRUE
    
    optimx.res <- do.call(optimx,
                          c(list(par = start, fn = objective, gr = gradient, control = control),
                            dots)
    )
  }
  
  
  index.cv <- which(optimx.res$convcode == 0)
  if(length(index.cv) == 0){
    res <- list(par = rep(NA, length(start)),
                objective = NA,
                convergence = 1,
                iterations = NA,
                evaluations = c("function" = NA, "gradient" = NA),
                message = "optimx has not converged \n"
    )
  }else if(length(index.cv) == 1){
    
    par_tempo <- unlist(lapply(1:length(start), function(x){optimx.res[[x]][index.cv]})) 
    names(par_tempo) <- names(start)
    
    res <- list(par = par_tempo,
                objective = optimx.res$value[index.cv],
                convergence = optimx.res$convcode[index.cv],
                iterations = optimx.res$niter[index.cv],
                evaluations = c("function" = optimx.res$fevals[index.cv], "gradient" = optimx.res$gevals[index.cv]),
                message = paste0("One algorithm has converged: ",attributes(optimx.res)$row.names[index.cv], " (value = ",optimx.res$value[index.cv],")\n")
    )
  }else{
    index.cvBest <- index.cv[which.max(optimx.res$value[index.cv])]
    par_tempo <- unlist(lapply(1:length(start), function(x){optimx.res[[x]][index.cvBest]})) 
    names(par_tempo) <- names(start)
    
    res <- list(par = par_tempo,
                objective = optimx.res$value[index.cvBest],
                convergence = optimx.res$convcode[index.cvBest],
                iterations = optimx.res$niter[index.cvBest],
                evaluations = c("function" = optimx.res$fevals[index.cvBest], "gradient" = optimx.res$gevals[index.cvBest]),
                message = paste0("Several algorithms have converged but only the best convergence point is retained \n",
                                 "method: ",paste(attributes(optimx.res)$row.names[index.cv], collapse = " | "), "\n",
                                 "value : ",paste(optimx.res$value[index.cv], collapse = " | "), "\n")
                
    )
  }
  
  return(res)
  
}