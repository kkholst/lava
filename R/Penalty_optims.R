proxGrad <- function(start, objective, gradient, hessian,...){
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain",
                  "proxGrad.method", "fast", "browser")
  dots <- list(...)
  control <- dots$control
  control <- control[names(control) %in% PGcontrols]
  
  penalty <- dots$penalty
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  n.coef <- length(start)
  
  index.penaltyCoef <- which(names(start) %in% penalty$names.penaltyCoef)
  
  #### definition of the operator
  if(all(penalty$lambda1 == 0) && all(penalty$lambda2 == 0)){
    proxOperator <- function(x, step, lambda1, lambda2){x}
  }else if(all(penalty$lambda2 == 0)){
    test.penalty <- 1:n.coef %in% index.penaltyCoef
    
    proxOperator <- function(x, step, lambda1, lambda2){
      mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty)
    }
  }else if(all(penalty$lambda1 == 0)){
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
       resLassoPath <- LassoPath_lvm(beta0 = start, objectiveLv = objective, gradientLv = gradient, hessianLv = hessian,
                                     indexPenalty = index.penaltyCoef, indexNuisance = which(names(start) %in% penalty$names.varCoef), 
                                     sd.X = penalty$sd.X, base.lambda1 = penalty$lambda1, lambda2 = penalty$lambda2, 
                                     fix.nuisance = dots$fix.sigma, proxOperator = proxOperator, control = control)
       
      names(resLassoPath) <- c("lambda1", "lambda2", names(start))
      
#     if(length(constrain)>0){
#       lambda2 <- rep(0, n.coef)
#       lambda1 <- rep(0, n.coef)
#       
#       for(iter_lambda in 1:nrow(resLassoPath)){
#         
#         iter_start  <- unlist(resLassoPath[iter_lambda,names(start)])
#         lambda1[index.penaltyCoef] <- resLassoPath[iter_lambda,"lambda1"]
#        
#         resTempo <- do.call(dots$proxGrad.method,
#                             list(start = iter_start, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
#                                  lambda1 = lambda1, lambda2 = penalty$lambda2, constrain = iter_start[names(iter_start) %in% names(constrain) == FALSE],
#                                  iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast))
#         
#         resLassoPath[iter_lambda,names(constrain)] <- resTempo$par[names(constrain)]
#       }
#     }
    if(control$constrain == TRUE){
      resLassoPath[,penalty$names.varCoef] <- exp(resLassoPath[,penalty$names.varCoef])
    }
    
    res <- list(par = start,
                convergence = 0,
                iterations = 0,
                evaluations = c("function" = 0, "gradient" = 0),
                message = resLassoPath,
                objective = NA)
  }else{
    
    ## one step lasso
    lambda1 <- rep(0, n.coef)
    lambda1[index.penaltyCoef] <- penalty$lambda1 
    lambda2 <- rep(0, n.coef)
    lambda2[index.penaltyCoef] <- penalty$lambda2 
 
    res <- do.call(dots$proxGrad.method,
                   list(start = start, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                        lambda1 = lambda1, lambda2 = lambda2, constrain = constrain,
                        iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))

    res$objective <- objective(res$par, penalty = penalty)
     
    if(!is.null(constrain)){
      res2 <- do.call(dots$proxGrad.method,
                      list(start = res$par, step = res$step, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                           lambda1 = lambda1, lambda2 = lambda2, constrain = res$par[names(res$par) %in% names(constrain) == FALSE],
                           iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast))
      res$par <- res2$par
    }
  }
  
  ### export
  return(res)
}


ISTA <- function(start, step = NULL, proxOperator, hessian, gradient, objective,
                 lambda1, lambda2, constrain, n.backtracking = 10, eta.backtracking = 0.5,
                 iter.max, abs.tol, rel.tol, fast, trace = FALSE){

  ## initialisation
  test.cv <- FALSE
  iter <- 1
  x_k <- start 
  if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  if(fast == TRUE){
    t_k <- 1
    y_k <- x_k
  }
  
  if(is.null(step)){
    step <- 1/max(abs(eigen(hessian(start, penalty = NULL))$value))
  }
  
  if(trace){
    .pb <- txtProgressBar(min = 0, max = iter.max, style = 3)
  }
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){
  
    if(trace){setTxtProgressBar(.pb, iter)}
    x_km1 <- x_k
    iter_back <- 0
    diff_back <- 1
    
    while( (iter_back < n.backtracking) && (diff_back > 0) ){
      stepTest <- step*eta.backtracking^iter_back
    
      if(fast == TRUE){
        t_km1 <- t_k
        x_k <- proxOperator(x = y_k - stepTest * gradient(y_k, penalty = NULL), 
                            step = stepTest, lambda1 = lambda1, lambda2 = lambda2)
        t_k <- (1 + sqrt(1 + 4 * t_km1^2)) / 2
        y_k <- x_k + (t_km1-1)/t_k * (x_k - x_km1)
      }else{
        x_k <- proxOperator(x = x_km1 - stepTest * gradient(x_km1, penalty = NULL), 
                            step = stepTest, lambda1 = lambda1, lambda2 = lambda2)
      }
      if(!is.null(constrain)){
        x_k[names(constrain)] <- constrain
      }
      
      diff_x <- x_k - x_km1
      current_obj <- objective(x_k)
      diff_obj <- current_obj - objective(x_km1)
      diff_back <- diff_obj - (crossprod(diff_x, gradient(x_km1)) + 1/(2*stepTest) * crossprod(diff_x))
      iter_back <- iter_back + 1
    }
     
    step <- 1/max(abs(eigen(hessian(x_k, penalty = NULL))$value))#stepTest
    iter <- iter + 1
    absDiff <- abs(diff_obj) < abs.tol  #(abs(diff_x) < abs.tol)
    relDiff <- abs(diff_obj)/abs(current_obj) < rel.tol #(abs(diff_x)/abs(x_k) < rel.tol)
    test.cv <- absDiff + relDiff > 0  #all(  absDiff + ifelse(is.na(relDiff),0,relDiff) > 0 )
    
    # cat(step," ",iter_back, " ", max(abs(diff_x))," ",diff_obj,"\n")
    
  }
  
  if(trace){cat("\n")}
  
  ## export
  message <- if(test.cv){"Sucessful convergence \n"
  }else{
    paste("max absolute/relative difference: ",max(absDiff),"/",max(relDiff)," for parameter ",which.max(absDiff),"/",which.max(relDiff),"\n")
  }
  
  return(list(par = x_k,
              step = step,
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