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
#' @param constrain parameters to be constrained at a given value
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
ISTA <- function(start, proxOperator, hessian, gradient, objective,
                 lambda1, lambda2, group.lambda1,
                 step, BT.n, BT.eta, constrain, 
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
    BT.n <- 1
    BT.eta <- 1
  }else{
    stepMin <- step*BT.eta^BT.n
  }
  
  obj_k <- try(objective(x_k))
  if("try-error" %in% class(obj_k)){obj_k <- Inf}
  
  test.cv <- FALSE
  iter <- 1
  
  ## loop
  while(test.cv == FALSE && iter <= iter.max){
    
    x_km1 <- x_k
    obj_km1 <- obj_k
    iter_back <- 0
    diff_back <- 1
    
    while( (iter_back < BT.n) && (diff_back > 0) ){
      stepBT <- step*BT.eta^iter_back
      
      if(fast == 1){ # FISTA
        grad_tempo <- gradient(y_k)
        
        x_k <- proxOperator(x = y_k - stepBT * grad_tempo, 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
        
        diff_tempo <- x_k - y_k
        obj_tempo <- objective(y_k)
        
      }else if(fast == 2){ # Monotone APG
        grad_tempo <- gradient(y_k)
        
        z_k <- proxOperator(x = y_k - stepBT * grad_tempo, 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        v_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        x_k <- if(objective(z_k)<=objective(v_k)){z_k}else{v_k}
        if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
        
        diff_tempo <- x_k - y_k
        grad_tempo <- gradient(y_k)
        obj_tempo <- objective(y_k)
        
      }else if(fast == 3){ # Nesterov step
        z_km1 <- z_k
        z_k <- proxOperator(x = x_km1 - stepBT * gradient(x_km1), 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        x_k <- z_km1 + (z_k - z_km1) * iter / (iter + 3)
        if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
        
        diff_tempo <- x_k - z_k#x_km1
        grad_tempo <- try(gradient(z_k))#(x_km1)#
        obj_tempo <- try(objective(z_k))#obj_km1#
        
      }else{ # ISTA
        grad_tempo <- gradient(x_km1)
       
        x_k <- proxOperator(x = x_km1 - stepBT * grad_tempo, 
                            step = stepBT, lambda1 = lambda1, lambda2 = lambda2, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2)
        
        if(!is.null(constrain)){x_k[names(constrain)] <- constrain}
        diff_tempo <- x_k - x_km1
        obj_tempo <- obj_km1
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
      
      step <- max(stepMin,min(stepMax, stepBT/sqrt(BT.eta)))
    }
    # could also be computed using the Barzilai-Borwein Method     
    # step = crossprod(diff_x) / crossprod(diff_x, gradient(x_km1) - gradient(x_k)) 
    
    iter <- iter + 1
    absDiff <- abs(obj_k - obj_km1) < abs.tol
    relDiff <- abs(obj_k - obj_km1)/abs(obj_k) < rel.tol
    test.cv <- absDiff + relDiff > 0
    
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