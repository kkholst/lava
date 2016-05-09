LassoPath_lvm <- function(beta0, objectiveLv, hessianLv, gradientLv, 
                          indexPenalty, indexNuisance,
                          sd.X, base.lambda1, lambda2, group.lambda1, 
                          step, n.BT, eta.BT, fix.nuisance, proxOperator, control,
                          gradientPen = NULL, iter_max = length(beta0)*2){
  
  p <- length(beta0)
  
  if(is.null(gradientPen)){
    gradientPen <- function(beta, grad_Lv, ...){
      ifelse(beta == 0, -sign(grad_Lv), -sign(beta))  # because grad_Lv in lava is -grad(Lv)
    }
  }
  
  #### initialisation
  M.beta <- matrix(0, nrow = 1, ncol = p)
  colnames(M.beta) <- names(beta0)
  M.beta[1,] <- beta0
  if(fix.nuisance == TRUE){
    if(control$constrain){
      M.beta[1,indexNuisance] <- 0
    }else{
      M.beta[1,indexNuisance] <- 1  
    }
  }
  
  V.lambda <- max( abs(-gradientLv(M.beta[1,]) * sd.X)[indexPenalty] )
  
  cv <- FALSE
  iter <- 1
  
  #### main loop
  while(iter < iter_max && cv == FALSE){
    cat("*")
    
    ## prediction step: next breakpoint assuming linear path and constant sigma
    resNode <- nextNode_lvm(hessianLv = hessianLv, gradientLv = gradientLv,  gradientPen = gradientPen,
                            beta = M.beta[iter,], lambda1 = V.lambda[iter] * base.lambda1,
                            indexPenalty = indexPenalty, indexNuisance = indexNuisance)
    
    newLambda <- V.lambda[iter] * (1 - resNode$gamma)
    
    if(newLambda < 0 || is.infinite(newLambda)){
      cv <- TRUE
    }else{
      # cat(newLambda," : ",paste(resNode$beta, collapse = " - "),"\n")
      
      if(any(lambda2>0)){ ## correction of the beta due to the L2 correction
        resNode$beta <- do.call("ISTA",
                                list(start = resNode$beta, proxOperator = proxOperator, hessian = hessianLv, gradient = gradientLv, objective = objectiveLv,
                                     lambda1 = newLambda*base.lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, 
                                     step = step, n.BT = n.BT, eta.BT = eta.BT, constrain = resNode$beta[indexNuisance],
                                     iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))$par
      }
      
      if(fix.nuisance == FALSE){ ## estimation of sigma
        warperObj <- function(sigma2){
          pp <- resNode$beta
          pp[indexNuisance] <- sigma2
          return(objectiveLv(pp))
        }
        
        warperGrad <- function(sigma2){
          pp <- resNode$beta
          pp[indexNuisance] <- sigma2
          return(sum(abs(gradientLv(pp)[indexNuisance])))
        }
        
        #       test1 <- optim(par = resNode$beta[indexNuisance], 
        #                      fn = warperObj,
        #                      gr = warperGrad,
        #                      method = "Brent", lower = 0, upper = M.beta[1,indexNuisance])$par
        
        suppressWarnings(
          resNode$beta[indexNuisance] <- optim(par = resNode$beta[indexNuisance], 
                                               fn = warperGrad, 
                                               lower = rep(1e-12,length(indexNuisance)), upper = M.beta[1, indexNuisance])$par
        )
        
        ## update lambda given that the product lambda*sigma must be a constant
        newLambda <- median(newLambda * M.beta[iter,names(resNode$beta[indexNuisance])] / resNode$beta[indexNuisance])
      }
      
      ## correction step
      resNode$beta <- do.call("ISTA",
                              list(start = resNode$beta, proxOperator = proxOperator, hessian = hessianLv, gradient = gradientLv, objective = objectiveLv,
                                   lambda1 = newLambda*base.lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, 
                                   step = step, n.BT = n.BT, eta.BT = eta.BT, constrain = if(fix.nuisance){resNode$beta[indexNuisance]}else{NULL},
                                   iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))$par
      
      V.lambda <- c(V.lambda, newLambda)
      M.beta <- rbind(M.beta, resNode$beta)
    } 
    
    iter <- iter + 1
  }
  cat("\n")
  
  #### export
  return(data.frame(lambda1 = V.lambda, lambda2 = NA, M.beta))
}

#' @title  Find the next value of the regularization parameter

nextNode_lvm <- function(hessianLv, gradientLv, gradientPen,
                         beta, lambda1, indexPenalty, indexNuisance){
  
  ##
  hess_Lv <- -hessianLv(beta)  # because hess_Lv in lava is -hess(Lv)
  grad_Lv <- -attr(hess_Lv, "grad")  # because grad_Lv in lava is -grad(Lv) // gradientLv(beta)
  grad_Pen <- gradientPen(beta, grad_Lv = grad_Lv)
  
  set_A <- union(setdiff(which(beta!=0),indexNuisance),
                 setdiff(which(abs(grad_Lv)/lambda1 > 1-1e-04),indexNuisance)
  )
  
  grad_B.A <- solve(hess_Lv[set_A,set_A, drop = FALSE], grad_Pen[set_A, drop = FALSE] * lambda1[set_A, drop = FALSE])
  grad_Rho <- hess_Lv[,set_A, drop = FALSE] %*% grad_B.A / lambda1 
  
  ## find next knot
  gamma1 <- -beta[set_A]/(grad_B.A*lambda1[set_A]) 
  gamma2Moins <- (1 - grad_Lv/lambda1)/(1 + grad_Rho)
  gamma2Moins[c(set_A,indexNuisance)] <- Inf
  gamma2Plus <- (1 + grad_Lv/lambda1)/(1 - grad_Rho)
  gamma2Plus[c(set_A,indexNuisance)] <- Inf
  
  # cat(paste(gamma1, collpase = " ")," | ",paste(gamma2Moins, collpase = " ")," | ",paste(gamma2Plus, collpase = " ")," | \n \n")
  
  test.min <- any(c(gamma1,gamma2Moins,gamma2Plus)>=0)
  if(test.min){
    gamma <- max(1e-04, min(c(gamma1,gamma2Moins,gamma2Plus)[c(gamma1,gamma2Moins,gamma2Plus)>=0]))
  }else{
    gamma <- Inf
  }
  
  ## update coef
  beta[set_A] <- beta[set_A] + gamma * grad_B.A
  
  return(list(gamma = gamma, beta = beta))
}


EPSODE <- function(beta, objective, gradient, hessian, V, indexPenalty, sweep = FALSE, 
                   step_lambda1, resolution_lambda1 = 1000, nstep_max = length(beta)*5, tol = 1e-7,
                   lambda2, group.lambda1, step, n.BT, eta.BT, proxOperator, control){
  
  #### preparation
  n.coef <- length(beta)
  lambda2_save <- lambda2
  lambda2 <- rep(0, n.coef)
  lambda2[indexPenalty] <- lambda2_save 
  envir <- environment()
  
  #### initialization
  iter <- 1
  if(step_lambda1 > 0){
    seq_lambda1 <- 0  
  }else{
    seq_lambda1 <-  max( abs(-gradient(beta) )[indexPenalty] )  
    beta <- do.call("ISTA",
                    list(start = beta, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                         lambda1 = seq_lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, 
                         step = step, n.BT = n.BT, eta.BT = eta.BT, constrain = NULL,
                         iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))$par
  }
  
  M.beta <- rbind(beta)
  setNE <- intersect(which(V %*% beta < 0),  indexPenalty)
  setZE <- intersect(which(V %*% beta == 0), indexPenalty)
  setPE <- intersect(which(V %*% beta > 0), indexPenalty)
  test.ncv <- TRUE
  
  #### main loop
  while(iter < nstep_max && test.ncv){
    
    ## current parameters
    iterLambda1 <- seq_lambda1[iter]
    iterBeta <- M.beta[iter,]
    if(length(setZE)>0){
      iterBeta[setZE] <- 0
    }
    
    cat(iter," ",iterLambda1," ", paste(iterBeta, collapse = " "),"\n")  
  
    ## Solve ODE 
    cat("A ")
    lambda.ode <- seq(iterLambda1, max(0, iterLambda1 + step_lambda1), length.out = resolution_lambda1)
    res.ode <- ode(y = c(0,iterBeta), 
                   times = unique(round(lambda.ode, digit = 6)), 
                   func = EPSODE_odeBeta,  
                   parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, indexPenalty = indexPenalty, 
                               step_lambda1 = step_lambda1, tol = tol,
                               sweep = sweep, envir = envir)
    )
    
    ## detect breakpoints
    res.breakpoint <- EPSODE_breakpoint(res.ode = res.ode, V = V, hessian = hessian, gradient = gradient,
                                        setNE = setNE, setZE = setZE, setPE = setPE, indexPenalty = indexPenalty)
    
    ## more precise estimate
    if(any(res.breakpoint$breakpoint.sign, res.breakpoint$breakpoint.constrain)){
      cat("B ")
      lambda.ode <- seq(res.ode[res.breakpoint$index_minus, 1], res.ode[res.breakpoint$index_plus, 1], length.out = resolution_lambda1)
      res.ode <- ode(y = c(0,res.ode[res.breakpoint$index_minus,-(1:2)]), 
                     times = unique(round(lambda.ode, digit = 6)), 
                     func = EPSODE_odeBeta,  
                     parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, indexPenalty = indexPenalty, 
                                 step_lambda1 = step_lambda1, tol = tol,
                                 sweep = sweep, envir = envir)
      )

      res.breakpoint <- EPSODE_breakpoint(res.ode = res.ode, V = V, hessian = hessian, gradient = gradient,
                                          setNE = setNE, setZE = setZE, setPE = setPE, indexPenalty = indexPenalty)

    }
    
    ## prepare update
    newLambda1 <- res.breakpoint$newLambda
    newBeta <- res.ode[res.breakpoint$newLambda,-(1:2)]
    
    ## correction step
    cat("C ")
    lambda1 <- rep(0, n.coef)
    lambda1[indexPenalty] <- newLambda1
    newBeta <- do.call("ISTA",
                       list(start = newBeta, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                            lambda1 = lambda1, lambda2 = lambda2, group.lambda1 = group.lambda1, 
                            step = step, n.BT = n.BT, eta.BT = eta.BT, constrain = NULL,
                            iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$fast, trace = FALSE))$par
    cat("\n ")
    setNE <- intersect(which(V %*% newBeta < -tol),  indexPenalty)
    setPE <- intersect(which(V %*% newBeta > tol), indexPenalty)
    setZE <- setdiff(indexPenalty, c(setNE,setPE))
    
    ## updates
    M.beta <- rbind(M.beta, newBeta)
    seq_lambda1 <- c(seq_lambda1, newLambda1)
    iter <- iter + 1
    if(step_lambda1 > 0){
      test.ncv <- (length(setNE) > 0 || length(setPE) > 0 )
    }else{
      test.ncv <- length(setZE)
    }
    
  }
  
  #### export
  return(data.frame(lambda1 = unname(seq_lambda1), lambda2 = NA, M.beta))
}

EPSODE_odeBeta <- function(t, y, ls.args){
  
 y <- y[-1]
  
  #### test knot
  if(length(ls.args$setNE)>0 && any(y[ls.args$setNE]>0)){ ## any negative parameter that becomes positive: stop algorithm
    index <- ls.args$setNE[which.max(y[ls.args$setNE])]
    return(list(c(1,rep(0,length(y)))))
  }
  if(length(ls.args$setPE)>0 && any(y[ls.args$setPE]<0)){ ## any positive parameter that becomes negative: stop algorithm
    index <- ls.args$setNE[which.min(y[ls.args$setNE])]
    return(list(c(1,rep(0,length(y)))))
  }
  
  #### estimate uz
  uz <- rep(0, length(y))
  if(length(ls.args$setNE)>0){
    uz <- uz - colSums(ls.args$V[ls.args$setNE,,drop = FALSE])
  }
  if(length(ls.args$setPE)>0){
    uz <- uz  + colSums(ls.args$V[ls.args$setPE,,drop = FALSE])
  }
  
  #### estimate Q and P
  if(ls.args$sweep == FALSE || t == 0){
    
    #     H <- ls.args$hessian(y) #### problem when using LAVA derivatives
    H <- -hessianO(y) ; attr(H, "grad") <- -gradientO(y) #
    #     if(min(eigen(H)$values) < 1e-5){browser()} #### badly conditionned hessian
    H_m1 <- solve(H)
    
    if(length(ls.args$setZE) == 0){
      R <- NULL
      Q <- NULL
      P <- H_m1
    }else{
      Uz <- ls.args$V[ls.args$setZE,,drop = FALSE]
      
      R <- solve(Uz %*% H_m1 %*% t(Uz))
      Q <- H_m1 %*% t(Uz) %*% R
      P <- H_m1 - Q %*% Uz %*% H_m1 
    }
    
  } else {
    
    #     res.ode <- ode(y = get(x = "P_hist", envir = ls.args$envir), 
    #                    times = c(get(x = "t_hist", envir = ls.args$envir),t), 
    #                    func = EPSODE_odeP,  
    #                    parm = list(hessian = ls.args$hessian, uz = uz, beta = y)
    #     )
    ### same for R and Q
    
  }
  
  ## check constrains
  if(!is.null(Q)){
    
    G <- attr(H, "grad") # -gradient(iterBeta)
    s <- - t(Q) %*% (1 / t * G + uz)
   
    if(any( abs(s) > 1)){
     
      # cat(t, " : ",paste(y, collapse = " "))
#       cat(t, " ")
#       print(H)
#       print(s)
    
      index <- ls.args$setZE[which.max(abs(s))]
      return(list(c(1,rep(0,length(y)))))
    }
    
  }
  
  
  
  #   cat(t, " : ", paste(y, collapse = " "),"\n")
  #   
  #   if(t > 0){
  #     get("H_lava", envir = ls.args$envir)
  #     get("H_m1_lava", envir = ls.args$envir)  
  #     get("P_lava", envir = ls.args$envir)  
  #     get("uz_lava", envir = ls.args$envir)  
  #     get("y_lava", envir = ls.args$envir)  
  #     get("t_lava", envir = ls.args$envir)  
  #     
  #     H_lava <<- c(H_lava, list(H))
  #     H_m1_lava <<- c(H_m1_lava, list(H_m1))
  #     P_lava <<- c(P_lava, list(P))
  #     uz_lava <<- c(uz_lava, list(uz))
  #     y_lava <<- c(y_lava,list(y))
  #     t_lava <<- c(t_lava,list(t))
  #     
  #   }else{
  #     H_lava <<- list(H)
  #     H_m1_lava <<- list(H_m1)
  #     P_lava <<- list(P)
  #     uz_lava <<- list(uz)
  #     y_lava <<- list(y)
  #     t_lava <<- list(t)
  #   }
  
  
  #   print(P)
  #   print(uz)
  #   browser()
  #   cat("*")
  ## export
  Puz <- P %*% uz
  Puz[ls.args$setZE] <- 0
  return(list(c(0,-Puz * sign(ls.args$step_lambda1))))
}

EPSODE_odeP <- function(t, y, ls.args){
  dH <- genD(ls.args$hessian, as.numeric(ls.args$beta))
  return( (y %x% y) %*% dH %*% y %*% ls.args$uz )
}


EPSODE_breakpoint <- function(res.ode, V, hessian, gradient,
                              setNE, setZE, setPE, indexPenalty){
  
  n.lambda <- nrow(res.ode)
  
  #### detect breakpoint
  index.breakpoint <- which(abs(res.ode[,2])>0) #which(duplicated(res.ode[,-(1:2)]))
  
  if(length(index.breakpoint) == 0){ ## no breakpoint
    res <- list(breakpoint.sign = FALSE,
                breakpoint.constrain = FALSE,
                index_ode = NULL,
                index_coef = NULL,
                newLambda = res.ode[n.lambda,1])
  }else{ 
    index.breakpoint <- index.breakpoint[1]
    index.changeSign <- apply(res.ode[,2+c(setNE,setPE), drop = FALSE], 2, function(x){which(diff(sign(x))!=0)[1]})
    
   if(!all(is.na(index.changeSign))){ ## change sign
     index.breakpoint <- min(index.changeSign, na.rm = TRUE) ## to be sure
     
     res <- list(breakpoint.sign = TRUE,
                 breakpoint.constrain = FALSE,
                 index_ode = index.breakpoint,
                 index_coef = which.min(index.changeSign),
                 newLambda = res.ode[index.breakpoint,1])
    }else{ ## no zero parameter
      y <-  res.ode[index.breakpoint,-(1:2)]
      
      uz <- rep(0, length(y))
      if(length(setNE)>0){
        uz <- uz - colSums(V[setNE,,drop = FALSE])
      }
      if(length(setPE)>0){
        uz <- uz  + colSums(V[setPE,,drop = FALSE])
      }
      #     H <- ls.args$hessian(y) #### problem when using LAVA derivatives
      H <- -hessianO(y) ; attr(H, "grad") <- -gradientO(y) #
      H_m1 <- solve(H)
      Uz <- V[setZE,,drop = FALSE]
      R <- solve(Uz %*% H_m1 %*% t(Uz))
      Q <- H_m1 %*% t(Uz) %*% R
      
      G <- attr(H, "grad") # -gradient(iterBeta)
      s <- - t(Q) %*% (1 / res.ode[index.breakpoint,1] * G + uz)
      index <- which.max(abs(s))
      
      res <- list(breakpoint.sign = FALSE,
                  breakpoint.constrain = TRUE,
                  index_ode = index.breakpoint,
                  index_coef = index,
                  newLambda = res.ode[index.breakpoint,1])
    }
  }
  
  #### post treatment
  res$index_plus <- min(n.lambda, res$index_ode + 1)
  res$index_minus <- max(1, res$index_ode - 1)
  
  #### export
  return(res)
  
}





