proxGrad <- function(start, objective, gradient, hessian,...){
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain",
                  "group.penaltyCoef",
                  "fast", "step", "n.BT", "eta.BT", "trace")
  dots <- list(...)
  control <- dots$control
  control <- control[names(control) %in% PGcontrols]
  
  penalty <- dots$penalty
  if (length(dots$trace)>0 && dots$trace>0) cat("\n")
  n.coef <- length(start)
  
  #### definition of the operator
  if(any(penalty$group.penaltyCoef>1)){ # grouped lasso
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){
      levels.penalty <- unique(test.penalty1)
      for(iter_group in 1:length(levels.penalty)){
        index_group <- which(test.penalty1==levels.penalty[iter_group])
        x[index_group] <- proxE2(x = x[index_group], step = step, lambda = lambda1[index_group])
      }
      return(x)
    }
  }else if(all(penalty$lambda1 == 0) && all(penalty$lambda2 == 0)){
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){x}
  }else if(all(penalty$lambda2 == 0)){
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){
      mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty1)
      #proxL1_cpp(x, step, lambda1, which(test.penalty1>0)-1, sum(test.penalty1>0))
    }
  }else if(all(penalty$lambda1 == 0)){
     proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){
      mapply(proxL2, x = x, step = step, lambda = lambda2, test.penalty = test.penalty2)
      #proxL2_cpp(x, step, lambda1, which(test.penalty2>0)-1, sum(test.penalty2>0))
    }
  }else{
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){
      mapply(proxL2, 
             x = mapply(proxL1, x = x, step = step,  lambda = lambda1, test.penalty = test.penalty1),
             step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
  }
  
  #### Penalty
  
  ## check
  if(length(setdiff(penalty$names.penaltyCoef, names(start)))>0){
    message("proxGrad: some penalty will not be applied because the corresponding parameter is used as a reference \n",
            "non-applied penalty: ",paste(setdiff(penalty$names.penaltyCoef, names(start)), collapse = " "),"\n")
    penalty$group.penaltyCoef <- penalty$group.penaltyCoef[penalty$names.penaltyCoef %in% names(start)]
    penalty$names.penaltyCoef <- penalty$names.penaltyCoef[penalty$names.penaltyCoef %in% names(start)]
  }
  index.penaltyCoef <- which(names(start) %in% penalty$names.penaltyCoef)

  ## grouped lasso
  tempo <- penalty$group.penaltyCoef
  penalty$group.penaltyCoef <- setNames(rep(0, n.coef), names(start))
  penalty$group.penaltyCoef[penalty$names.penaltyCoef ] <- tempo
  
  ## TO BE REMOVED
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
  
#     V <- matrix(0, nrow = length(start), ncol = length(start))
#     diag(V)[index.penaltyCoef] <- 1
#     EPSODE(gradient = gradient, hessian = hessian, beta = start, V = V, lambda = penalty$lambda1, indexPenalty = index.penaltyCoef)
    
     ## path regularization
       resLassoPath <- LassoPath_lvm(beta0 = start, objectiveLv = objective, gradientLv = gradient, hessianLv = hessian,
                                     indexPenalty = index.penaltyCoef, indexNuisance = which(names(start) %in% penalty$names.varCoef), 
                                     sd.X = penalty$sd.X, base.lambda1 = penalty$lambda1, lambda2 = penalty$lambda2, group.lambda1 = penalty$group.penaltyCoef,
                                     step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT, 
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
   
    res <- ISTA(start = start, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                lambda1 = lambda1, lambda2 = lambda2, group.lambda1 = penalty$group.penaltyCoef, constrain = constrain,
                step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT, trace = if(!is.null(control$trace)){control$trace}else{FALSE},
                iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast)

    res$objective <- objective(res$par, penalty = penalty)
     
    if(!is.null(constrain)){
      res2 <- ISTA(start = res$par, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                   lambda1 = lambda1, lambda2 = lambda2, group.lambda1 = penalty$group.penaltyCoef, constrain = res$par[names(res$par) %in% names(constrain) == FALSE],
                   step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT,
                   iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast)
      res$par <- res2$par
    }
  }
  
  ### export
  return(res)
}


ISTA <- function(start, proxOperator, hessian, gradient, objective,
                 lambda1, lambda2, group.lambda1,
                 step, n.BT, eta.BT,  constrain, 
                 iter.max, abs.tol, rel.tol, fast, trace = FALSE){

   ## initialisation
  test.penalty1 <- group.lambda1
  test.penalty2 <- (lambda2>0)
    
  x_k <- start 
  if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
  if(fast %in% 1:2){
    t_k <- 1
    y_k <- x_k
  }else if(fast == 3){
    z_k <- x_k
  }
   
    stepMax <- step 
  if(is.null(step)){
    step <- 1/max(abs(eigen(hessian(start, penalty = NULL))$value))
    n.BT <- 1
    eta.BT <- 1
  }else{
    if(is.null(n.BT) || is.null(eta.BT) ){
      stop("ISTA: if argument \'step\' is not null then argments \'n.BT\' and \'eta.BT\' must be specified \n")
    }
     stepMin <- step*eta.BT^n.BT
    # step <- 1/max(abs(eigen(hessian(start, penalty = NULL))$value))
  }
  
  obj_k <- try(objective(x_k))
  if("try-error" %in% class(obj_k)){obj_k <- Inf}
  
  test.cv <- FALSE
  iter <- 1
  
#   if(trace){
#     .pb <- txtProgressBar(min = 0, max = iter.max, style = 3)
#   }
  
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){

    # if(trace){setTxtProgressBar(.pb, iter)}
    x_km1 <- x_k
    obj_km1 <- obj_k
    iter_back <- 0
    diff_back <- 1
    
     while( (iter_back < n.BT) && (diff_back > 0) ){
      stepBT <- step*eta.BT^iter_back
      
      if(fast == 1){ # Bech and Teboulle (2009 A Fast Iterative Shrinkage-Thresholding Algorithm)
        x_k <- proxOperator(x = y_k - stepBT * gradient(y_k, penalty = NULL), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
         
        diff_tempo <- x_k - y_k
        grad_tempo <- gradient(y_k)
        obj_tempo <- objective(y_k)
        
      }else if(fast == 2){ # Li Accelerated Proximal Gradient Methods
        z_k <- proxOperator(x = y_k - stepBT * gradient(y_k, penalty = NULL), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        v_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1, penalty = NULL), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        x_k <- if(objective(z_k)<=objective(v_k)){z_k}else{v_k} # fast <- FALSE is not in the original method
        
        diff_tempo <- x_k - y_k
        grad_tempo <- gradient(y_k)
        obj_tempo <- objective(y_k)
        
      }else if(fast == 3){ # Nesterov step A SPARSE-GROUP LASSO
        z_km1 <- z_k
        z_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1, penalty = NULL), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        x_k <- z_km1 + (z_k - z_km1) * iter / (iter + 3)
        diff_tempo <- x_k - z_k#x_km1
        grad_tempo <- try(gradient(z_k))#(x_km1)#
        obj_tempo <- try(objective(z_k))#obj_km1#
        
      }else{ # normal step
        x_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1, penalty = NULL), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        diff_tempo <- x_k - x_km1
        grad_tempo <- gradient(x_km1)
        obj_tempo <- obj_km1
      }
      
      
      if(!is.null(constrain)){
        x_k[names(constrain)] <- constrain
      }
      
      obj_k <- try(objective(x_k))
      if(any("try-error" %in% c(class(obj_k), class(obj_tempo), class(grad_tempo)))){
        diff_back <- Inf
      }else{
        # crossprod is not supposed to be negative ...
        diff_back <- obj_k - (obj_tempo + crossprod(diff_tempo, grad_tempo) + 1/(2*stepBT) * crossprod(diff_tempo) ) 
      }
      iter_back <- iter_back + 1
      
    }
    
    if(fast == 1){
      t_km1 <- t_k
      t_k <- (1 + sqrt(1 + 4 * t_km1^2)) / 2
      y_k <- x_k + (t_km1-1)/t_k * (x_k - x_km1)
    }else if(fast == 2){
      t_km1 <- t_k
      t_k <- (1 + sqrt(1 + 4 * t_km1^2)) / 2
      y_k <- x_k + t_km1/t_k * (z_k - x_k) + (t_km1-1)/t_k * (x_k - x_km1)
    }
    
    if(is.null(stepMax)){
      step <- 1/max(abs(eigen(hessian(x_k, penalty = NULL))$value)) 
    }else{
      
      step <- max(stepMin,min(stepMax, stepBT/sqrt(eta.BT)))
    }
    # could also be computed using the Barzilai-Borwein Method     
    # step = crossprod(diff_x) / crossprod(diff_x, gradient(x_km1) - gradient(x_k)) 
   
    iter <- iter + 1
    absDiff <- abs(obj_k - obj_km1) < abs.tol  #(abs(diff_x) < abs.tol)
    relDiff <- abs(obj_k - obj_km1)/abs(obj_k) < rel.tol #(abs(diff_x)/abs(x_k) < rel.tol)
    test.cv <- absDiff + relDiff > 0  #all(  absDiff + ifelse(is.na(relDiff),0,relDiff) > 0 )
     
    if(trace){cat(stepBT," ",iter_back, " ", max(abs(x_k - x_km1))," ",obj_k - obj_km1,"\n")}
  }
  
  if(trace){cat("\n")}
  
  ## export
  message <- if(test.cv){"Sucessful convergence \n"
  }else{
    paste("max absolute/relative difference: ",max(absDiff),"/",max(relDiff)," for parameter ",which.max(absDiff),"/",which.max(relDiff),"\n")
  }
  
  return(list(par = x_k,
              step = stepBT,
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