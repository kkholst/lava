ISTA <- function(start, objective, gradient, hessian,...){
  
  ISTAcontrols <- c("iter.max","trace","abs.tol","rel.tol","FAST")
  dots <- list(...)
  control <- dots$control
  control <- control[names(control) %in% ISTAcontrols]
  penalty <- dots$penalty
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  
  ## definition of the operators
  n.coef <- length(start)
  
  proxL1 <- function(x, lambda, test.penalty){
    if(test.penalty){
      max(0, (abs(x) - lambda*t_k)) * sign(x)
    }else{
      x
    }
  }
  
  #### initialisation
  test.cv <- FALSE
  iter <- 1
  x_k <- start 
  
  t_k <- 1/max(abs(eigen(hessian(x_k, penalty = penalty))$value)) # may be computationnaly expensive 
  lambda1 <- rep(0, length(x_k))
  lambda1[penalty$index.coef] <- penalty$lambda1 
  lambda2 <- rep(0, length(x_k))
  lambda2[penalty$index.coef] <- penalty$lambda2 
  
  ## requested operator
  if(all(lambda1 == 0)){
    proxOperator <- function(x, lambda1){x}
  }else{
    test.penalty <- lambda1>0
    
    proxOperator <- function(x, lambda1){
      mapply(proxL1, x = x, lambda = lambda1, test.penalty = test.penalty)
    }
  }
 
  #### main
  if(("FAST" %in% names(control) == FALSE) || (control$FAST == FALSE)){
    ## ISTA
  
    while(test.cv == FALSE && iter <= control$iter.max){
      x_km1 <- x_k
      x_k <- proxOperator(x = x_km1 - t_k * gradient(x_km1, penalty = penalty), 
                          lambda1 = lambda1)
      
      iter <- iter + 1
      absDiff <-  (abs(x_k-x_km1) < control$abs.tol)
      relDiff <- (abs(x_k-x_km1)/abs(x_k) < control$rel.tol)
      test.cv <- all(  absDiff + ifelse(is.na(relDiff),0,relDiff) > 0 )
     
    }
    
  }else{
    
    ## FISTA
    y_k <- x_k  
    t2_k <- 1
    
    while(test.cv==FALSE && iter <= control$iter.max){
      
      x_km1 <- x_k
      t2_km1 <- t2_k
      
      x_k <- proxOperator(x = y_k - t_k * gradient(y_k, penalty = penalty), 
                          lambda1 = lambda1)
      
      t2_k <- (1 + sqrt(1 + 4 * t2_k^2)) / 2
      y_k <- x_k + (t2_km1-1)/t2_k * (x_k - x_km1)
      
      iter <- iter + 1
      absDiff <-  (abs(x_k-x_km1) < control$abs.tol)
      relDiff <- (abs(x_k-x_km1)/abs(x_k) < control$rel.tol)
      test.cv <- all(  absDiff + ifelse(is.na(relDiff),0,relDiff) > 0 )
      
    }
  }
  
  ### export
  res <- list(par = x_k,
              objective = objective(x_k, penalty = penalty),
              convergence = as.numeric(test.cv==FALSE),
              iterations = iter,
              evaluations = c("function" = 0, "gradient" = iter)
              )      
  return(res)
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