#' @title step of a proximal gradient algorithm
#' @param start initial values for the parameters
#' @param proxOperator proximal operator corresponding to the penalization applied to the log likelihood
#' @param hessian second derivative of the likelihood given by lava. Only used to estimate the step parameter of the algorithm when step = NULL
#' @param gradient first derivative of the likelihood given by lava. 
#' @param objective likelihood given by lava. Used to adjust the step parameter when using backtracking
#' @param lambda1 L1 penalization parameter
#' @param lambda2 L2 penalization parameter
#' @param group.lambda1 group to which each parameter belongs. 0 mean individual lasso otherwise parameters are groupes according to their group.lambda1 value
#' @param step maximun step for the proximal gradient algorithm. 
#' If NULL the step is estimated using the inverse of the maximal eigenvalue of the hessian (in absolute value) and re-estimated at each step
#' Otherwise backtracking is used.
#' @param BT.n number of backtracking steps
#' @param BT.eta multiplicative factor for the step 
#' @param iter.max maximum number of iterations
#' @param abs.tol convergence is the difference in likelihood between two consecutive steps is below this threshold
#' @param rel.tol convergence is the relative difference  in likelihood between two consecutive steps is below this threshold
#' @param fast type of iteration
#' 0 correspond to the ISTA step as described in Bech 2009
#' 1 correspond to the FISTA step as described in Bech 2009
#' 2 correspond to the Monotone APG as described in Li 2015
#' 3 correspond to the Nesterov step as described in Simon 2013
#' @param trace should the convergence diagnostics be displayed at each step
#' 
#' @references 
#' Bech and Teboulle - 2009 A Fast Iterative Shrinkage-Thresholding Algorithm
#' Li 2015 - Accelerated Proximal Gradient Methods for Nonconvex Programming
#' Simon 2013 - A sparse group Lasso
proxGrad <- function(start, proxOperator, hessian, gradient, objective,
                     step, BT.n, BT.eta, force.descent, 
                     iter.max, abs.tol, rel.tol, method, trace = FALSE){
  
#   if(is.null(step)){
#     step <- 100/sum( (gradient(start)/start)^2)#add variance and x
#   }
  stepMax <- step 
  stepMin <- step*BT.eta^BT.n
  fct_errorLv <- function(e){warning("unable to compute the value of the likelihood - return Inf \n");return(Inf)}
  ## initialisation
  x_k <- start 

  obj.x_k <- tryCatch(objective(x_k), error = fct_errorLv)
  if(is.na(obj.x_k)){obj.x_k <- Inf}
  grad.x_k <- try(gradient(x_k))
  
  t_k <- t_kp1 <- if(method %in% c("FISTA","mFISTA")){1}else{NA}
  y_k <- if(method %in% c("FISTA","mFISTA")){x_k}else{NA} 
  z_k <- if(method == "Nesterov"){z_k}else{NA}
  obj.y_k <- if(method %in% c("FISTA","mFISTA")){obj.x_k}else{NA} 
  grad.y_k <- if(method %in% c("FISTA","mFISTA")){grad.x_k}else{NA} 
  
  test.cv <- FALSE
  iter <- 0
  iterAll <- 0
  
  if(trace>0){cat("stepBT"," ","iter_back", " ", "max(abs(x_k - x_km1))"," ","obj.x_k - obj.x_km1","\n")}
  
    ## loop
  while(test.cv == FALSE && iter <= iter.max){
 
    iter_back <- 0
    diff_back <- 1
    obj.x_kp1 <- +Inf
    
    while( (iter_back < BT.n) && (is.infinite(obj.x_kp1) || diff_back > 0) ){ # obj.x_kp1 > obj.x_k should not be needed
      stepBT <- step*BT.eta^iter_back
      iterAll <- iterAll + 1
      
      if(method == "ISTA"){
        res <- ISTA(x_k = x_k, obj.x_k = obj.x_k, grad.x_k = grad.x_k, 
                    proxOperator = proxOperator, step = stepBT)
      }else if(method == "FISTA"){
        res <- FISTA(x_k = x_k, y_k = y_k, t_kp1 = t_kp1, t_k = t_k,
                     obj.y_k = obj.y_k, grad.y_k = grad.y_k,
                     proxOperator = proxOperator, step = stepBT)
      }else if(method == "mFISTA"){
        res <- mFISTA(x_k = x_k, y_k = y_k, t_kp1 = t_kp1, t_k = t_k, 
                      obj.x_k = obj.x_k, obj.y_k = obj.y_k, grad.x_k = grad.x_k, grad.y_k = grad.y_k,
                      proxOperator = proxOperator, step = stepBT, objective = objective)
      }
      
      obj.x_kp1 <- tryCatch(objective(res$x_kp1), error = fct_errorLv) 
      if(is.na(obj.x_kp1)){obj.x_kp1 <- Inf}
      
      #diff_back <- max(obj.x_kp1 - res$Q, obj.x_kp1 - obj.x_k)
      if(force.descent == TRUE){
        diff_back <- obj.x_kp1 - obj.x_k
      }else{
        diff_back <- obj.x_kp1 - res$Q
      }
      
      iter_back <- iter_back + 1
    
      # cat("obj.x_kp1:",obj.x_kp1," | obj.x_k:",obj.x_k, " | res$Q:",res$Q,"\n")
    }
    
    if(force.descent && obj.x_kp1 > obj.x_k){break}
    
    absDiff <- abs(obj.x_kp1 - obj.x_k) < abs.tol
    relDiff <- abs(obj.x_kp1 - obj.x_k)/abs(obj.x_kp1) < rel.tol
    test.cv <- (absDiff + relDiff > 0)
    
    #### update
    x_km1 <- x_k
    obj.x_km1 <- obj.x_k
    grad.x_km1 <- grad.x_k
    
    x_k <- res$x_kp1
    obj.x_k <- obj.x_kp1
    grad.x_k <- try(gradient(res$x_kp1))
    
    if(method %in% c("FISTA","mFISTA")){
      y_km1 <- y_k
      y_k <- res$y_kp1
      t_k <- t_kp1
      
      t_kp1 <- (1 + sqrt(1 + 4 * t_k^2)) / 2
      obj.y_k <- tryCatch(objective(y_k), error = fct_errorLv)
      if(is.na(obj.y_k)){obj.y_k <- Inf}
      grad.y_k <- try(gradient(y_k))
    }
    step <- min(stepMax, stepBT/sqrt(BT.eta))
    
    iter <- iter + 1 
    
    if(trace>0){cat("|",stepBT," ",iter_back, " ", max(abs(x_k - x_km1))," ",obj.x_k - obj.x_km1,"\n")}
  }
  if(trace>0){cat("\n")}
  
  ## export
  message <- if(test.cv){"Sucessful convergence \n"
  }else{
    paste("max absolute/relative difference: ",max(abs(obj.x_k - obj.x_km1)),"/",max(abs(obj.x_k - obj.x_km1)/abs(obj.x_k))," for parameter ",which.max(absDiff),"/",which.max(relDiff),"\n")
  }
  
  return(list(par = x_k,
              step = stepBT,
              convergence = as.numeric(test.cv==FALSE),
              iterations = iter,
              iterationsAll = iterAll,
              evaluations = c("function" = 0, "gradient" = iter),
              message = message
  ))
}



#' @title Estimate an upper bound of obj.x
Qbound <- function(diff.xy, obj.y, grad.y, L){
  
  return(obj.y + crossprod(diff.xy, grad.y) + L/2 * crossprod(diff.xy))
  
}

ISTA <- function(x_k, obj.x_k, grad.x_k,
                 proxOperator, step){
  
  ## Step
  x_kp1 <- proxOperator(x = x_k - step * grad.x_k, step = step)
  
  ## Upper bound for backtracking
  Q <- Qbound(diff.xy = x_kp1 - x_k, obj.y = obj.x_k, grad.y = grad.x_k, L = 1/step)
  
  return(list(x_kp1 = x_kp1,
              Q = Q))
}

FISTA <- function(x_k, y_k, t_kp1, t_k, obj.y_k, grad.y_k,
                  proxOperator, step){
  
  ## step
  x_kp1 <- proxOperator(x = y_k - step * grad.y_k, step = step)
  y_kp1 <- x_kp1 + (t_k-1)/t_kp1 * (x_kp1 - x_k)
  
  ## Upper bound for backtracking
  Q <- Qbound(diff.xy = x_kp1 - y_k, obj.y = obj.y_k, grad.y = grad.y_k, L = 1/step)
  
  ## export
  return(list(x_kp1 = x_kp1,
              y_kp1 = y_kp1,
              Q = Q))
}

mFISTA <- function(x_k, t_kp1, t_k, y_k, obj.x_k, obj.y_k, grad.x_k, grad.y_k,
                   proxOperator, step, objective){
  
  ## step
  z_kp1 <- proxOperator(x = y_k - step * grad.y_k, step = step)
  Qz <- Qbound(diff.xy = z_kp1 - y_k, obj.y = obj.y_k, grad.y = grad.y_k, L = 1/step)
  
  resV <- ISTA(x_k = x_k, obj.x_k = obj.x_k, grad.x_k = grad.x_k,
               proxOperator = proxOperator, step = step)
  v_kp1 <- resV$x_kp1
  Qv <- resV$Q
  
  ## Upper bound for backtracking
  if(Qz <= Qv){
    x_kp1 <- z_kp1
    Q <- Qz
  }else{
    x_kp1 <- v_kp1
    Q <- Qv
  }
  
  ## extrapolation
  y_kp1 <- x_kp1 + t_k/t_kp1 * (z_kp1 - x_kp1) + (t_k-1)/t_kp1 * (x_kp1 - x_k)
  
  ## export
  return(list(x_kp1 = x_kp1,
              y_kp1 = y_kp1,
              Q = Q))
}

# if(fast == 2){ # Monotone APG
#   grad_tempo <- gradient(y_k)
#   
#    
#   x_k <- if(objective(z_k)<=objective(v_k)){z_k}else{v_k}
#   
#   diff_tempo <- x_k - x_km1#x_k - y_k
#   grad_tempo <- gradient(x_km1)#gradient(y_k)
#   obj_tempo <- obj_km1#objective(y_k)
#   
# }else if(fast == 3){ # Nesterov step
#   z_km1 <- z_k
#   z_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1), 
#                       step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
#   x_k <- z_km1 + (z_k - z_km1) * iter / (iter + 3)
#   
#   diff_tempo <- x_k - z_k#x_km1
#   grad_tempo <- try(gradient(z_k))#(x_km1)#
#   obj_tempo <- try(objective(z_k))#obj_km1#
#   
# }