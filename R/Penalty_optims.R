proxGrad <- function(start, objective, gradient, hessian,...){
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain",
                  "step_lambda1", "correction.step",
                  "fast", "step", "n.BT", "eta.BT", "trace")

  dots <- list(...)
  control <- dots$control
  control <- control[names(control) %in% PGcontrols]
  penalty <- dots$control$penalty
  dots$control$penalty <- NULL
  regularizationPath <- dots$control$regularizationPath
  dots$control$regularizationPath <- NULL
  fix.sigma <- dots$control$fix.sigma # to be removed
  dots$control$fix.sigma <- NULL # to be removed
  n.coef <- length(start)

  #### definition of the operator
  if(any(penalty$group.penaltyCoef>=1)){ # grouped lasso
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
  
  #### main
  if( regularizationPath > 0 ){
         
     ## path regularization
    if( regularizationPath == 1){
      
      resLassoPath <- LassoPath_lvm(beta0 = start, objectiveLv = objective, gradientLv = gradient, hessianLv = hessian,
                                    indexPenalty = index.penaltyCoef, indexNuisance = which(names(start) %in% penalty$names.varCoef), 
                                    sd.X = penalty$sd.X, base.lambda1 = penalty$lambda1, lambda2 = penalty$lambda2, group.lambda1 = penalty$group.penaltyCoef,
                                    step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT, 
                                    fix.nuisance = fix.sigma, proxOperator = proxOperator, control = control)
      
    }else{
      
      V <- matrix(0, nrow = length(start), ncol = length(start))
      diag(V)[index.penaltyCoef] <- 1
      
      resLassoPath <- EPSODE(beta = start, objective = objective, gradient = gradient, hessian = hessian, V = V, indexPenalty = index.penaltyCoef, 
                             step_lambda1 = control$step_lambda1, correction.step = control$correction.step,
                             lambda2 = penalty$lambda2, group.lambda1 = penalty$group.penaltyCoef,
                             step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT, proxOperator = proxOperator, 
                             control = control)
      
    }
    names(resLassoPath) <- c("lambda1", "lambda2", names(start))
    
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
    
    ## TO BE REMOVED
    if(fix.sigma){
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

    ## one step lasso
    lambda1 <- rep(0, n.coef)
    lambda1[index.penaltyCoef] <- penalty$lambda1 
    lambda2 <- rep(0, n.coef)
    lambda2[index.penaltyCoef] <- penalty$lambda2 
   
    res <- ISTA(start = start, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                lambda1 = lambda1, lambda2 = lambda2, group.lambda1 = penalty$group.penaltyCoef, constrain = constrain,
                step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT, trace = if(!is.null(control$trace)){control$trace}else{FALSE},
                iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast)
    
    res$objective <- objective(res$par) + penalized_objectivePen.lvm(res$par, 
                                                                     lambda1 = lambda1, lambda2 = lambda2)
     
    if(!is.null(constrain)){
      res2 <- ISTA(start = res$par, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                   lambda1 = lambda1, lambda2 = lambda2, group.lambda1 = penalty$group.penaltyCoef, constrain = res$par[names(res$par) %in% names(constrain) == FALSE],
                   step = control$step, n.BT = control$n.BT, eta.BT = control$eta.BT, trace = if(!is.null(control$trace)){control$trace}else{FALSE},
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
    step <- 1/max(abs(eigen(hessian(start))$value))
    n.BT <- 1
    eta.BT <- 1
  }else{
    if(is.null(n.BT) || is.null(eta.BT) ){
      stop("ISTA: if argument \'step\' is not null then argments \'n.BT\' and \'eta.BT\' must be specified \n")
    }
     stepMin <- step*eta.BT^n.BT
    # step <- 1/max(abs(eigen(hessian(start))$value))
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
        grad_tempo <- gradient(y_k)
        
        x_k <- proxOperator(x = y_k - stepBT * grad_tempo, 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
         
        diff_tempo <- x_k - y_k
        obj_tempo <- objective(y_k)
        
      }else if(fast == 2){ # Li Accelerated Proximal Gradient Methods
        grad_tempo <- gradient(y_k)
        
        z_k <- proxOperator(x = y_k - stepBT * grad_tempo, 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        v_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        x_k <- if(objective(z_k)<=objective(v_k)){z_k}else{v_k} # fast <- FALSE is not in the original method
        
        diff_tempo <- x_k - y_k
        grad_tempo <- gradient(y_k)
        obj_tempo <- objective(y_k)
        
      }else if(fast == 3){ # Nesterov step A SPARSE-GROUP LASSO
        z_km1 <- z_k
        z_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        x_k <- z_km1 + (z_k - z_km1) * iter / (iter + 3)
        diff_tempo <- x_k - z_k#x_km1
        grad_tempo <- try(gradient(z_k))#(x_km1)#
        obj_tempo <- try(objective(z_k))#obj_km1#
        
      }else{ # normal step
        grad_tempo <- gradient(x_km1)
      
        x_k <- proxOperator(x = x_km1 - stepBT * grad_tempo, 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        diff_tempo <- x_k - x_km1
        obj_tempo <- obj_km1
      }
      
      if(!is.null(constrain)){
        x_k[names(constrain)] <- constrain
      }
      
      obj_k <- try(objective(x_k))
      if(any("try-error" %in% c(class(obj_k), class(obj_tempo), class(grad_tempo)))){
        diff_back <- Inf
      }else{
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
      step <- 1/max(abs(eigen(hessian(x_k))$value)) 
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