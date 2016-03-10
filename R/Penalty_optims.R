proxGrad <- function(start, objective, gradient, hessian,...){
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol","proxGrad.method","fix.sigma","test.path")
  dots <- list(...)
  control <- dots$control
  control <- control[names(control) %in% PGcontrols]
  penalty <- dots$penalty
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  n.coef <- length(start)

  #### definition of the operator
  if(penalty$lambda1 == 0 && penalty$lambda2 == 0){
    proxOperator <- function(x, step, lambda1, lambda2){x}
  }else if(penalty$lambda2 == 0){
    test.penalty <- 1:n.coef %in% penalty$index.coef
    
    proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty)
    }
  }else if(penalty$lambda1 == 0){
    test.penalty <- 1:n.coef %in% penalty$index.coef
   
     proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL2, x = x, step = step, lambda = lambda2, test.penalty = test.penalty)
    }
  }else{
    test.penalty1 <- 1:n.coef %in% penalty$index.coef
    test.penalty2 <- 1:n.coef %in% penalty$index.coef
    
    proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL2, 
             x = mapply(proxL1, x = x, step = step,  lambda = lambda1, test.penalty = test.penalty1),
             step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
  }
  
  #### initialisation
  step <- 1/max(abs(eigen(hessian(start, penalty = NULL))$value)) # may be computationnaly expensive 
  
  lambda1 <- rep(0, n.coef)
  lambda1[penalty$index.coef] <- penalty$lambda1 
  lambda2 <- rep(0, n.coef)
  lambda2[penalty$index.coef] <- penalty$lambda2 
   
  #### main
  if(control$fix.sigma){
    names.constrain <- names(start)[grep(pattern = ",", names(start), fixed = TRUE)]
    constrain <- rep(1-dots$control$constrain, length(names.constrain))
    names(constrain) <- names.constrain
  }else{
    constrain <- NULL
  }
 
  res <- do.call(control$proxGrad.method,
                 list(start = start, step = step, proxOperator = proxOperator, gradient = gradient, 
                      lambda1 = lambda1, lambda2 = lambda2, constrain = constrain,
                      iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol))
  res$objective <- objective(res$par, penalty = penalty)
 
   #### path regularization
  if("test.path" %in% names(dots$control) && dots$control$test.path == TRUE){
    
    cat("current regularization: ",penalty$lambda1 ,"\n")
    
    pathRegularization(Y = as.matrix(df.data[,1,drop = FALSE]), 
                       X = cbind(1,as.matrix(df.data[,-1])),
                       start = start, objective = objective, gradient = gradient, hessian = hessian,
                       index.penalized = penalty$index.coef, constrain = constrain,
                       step = step, proxOperator = proxOperator, control = control, lambda2 = lambda2)
   
  }
   ### export

  return(res)
}



ISTA <- function(start, step = NULL, proxOperator, gradient, lambda1, lambda2, constrain,
                 iter.max, abs.tol, rel.tol){
  
  
  ## initialisation
  test.cv <- FALSE
  iter <- 1
  x_k <- start 
  if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){
    x_km1 <- x_k
    
    x_k <- proxOperator(x = x_km1 - step * gradient(x_km1, penalty = NULL), 
                        step = step, lambda1 = lambda1, lambda2 = lambda2)
    if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
    
    iter <- iter + 1
    absDiff <-  (abs(x_k-x_km1) < abs.tol)
    relDiff <- (abs(x_k-x_km1)/abs(x_k) < rel.tol)
    test.cv <- all(  absDiff + ifelse(is.na(relDiff),0,relDiff) > 0 )
    
  }
  
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

FISTA <- function(start, step = NULL, proxOperator, gradient, lambda1, lambda2, constrain,
                  iter.max, abs.tol, rel.tol){
  
  ## initialisation
  test.cv <- FALSE
  iter <- 1
  x_k <- start 
  y_k <- start  
  t_k <- 1
  if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){
    x_km1 <- x_k
    t_km1 <- t_k
   
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