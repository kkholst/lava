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
    # if(iter == 11){browser()}
    ## Solve ODE 
    cat("A ")

    lambda.ode <- seq(iterLambda1, max(0, iterLambda1 + step_lambda1), length.out = resolution_lambda1)
    res.ode <- ode(y = c(0,iterBeta), 
                   times = unique(round(lambda.ode, digit = 6)), 
                   func = EPSODE_odeBeta, method = "euler",
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
                     func = EPSODE_odeBeta, method = "euler",
                     parm = list(hessian = hessian, Vpen = V, setNE = setNE, setZE = setZE, setPE = setPE, indexPenalty = indexPenalty, 
                                 step_lambda1 = step_lambda1, tol = tol,
                                 sweep = sweep, envir = envir)
      )
  
      res.breakpoint <- EPSODE_breakpoint(res.ode = res.ode, V = V, hessian = hessian, gradient = gradient,
                                          setNE = setNE, setZE = setZE, setPE = setPE, indexPenalty = indexPenalty)
      
      if(res.breakpoint$breakpoint.sign == TRUE){
        setNE <- setdiff(setNE,  res.breakpoint$index_coef)
        setPE <- setdiff(setPE, res.breakpoint$index_coef)
        setZE <- union(setZE, res.breakpoint$index_coef)
      }
    
      if(res.breakpoint$breakpoint.constrain == TRUE){
        setZE <- setdiff(setZE, res.breakpoint$index_coef)
        if(res.breakpoint$sign_coef>0){
          setPE <- union(setPE, res.breakpoint$index_coef) 
        }else{
          setNE <- union(setNE,  res.breakpoint$index_coef)  
        }
      }
    }
    
    ## prepare update
    newLambda1 <- res.breakpoint$newLambda
    newBeta <- res.ode[res.breakpoint$index_ode,-(1:2)]
     
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
  return(data.frame(lambda1 = unname(seq_lambda1), lambda2 = NA, unname(M.beta)))
}

EPSODE_odeBeta <- function(t, y, ls.args){
  
  #### check cv
  if(abs(y[1])>0){
    return(list(rep(0,length(y))))
  }
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
    
    # H <- ls.args$hessian(y) #### problem when using LAVA derivatives
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
      # browser()
      index <- ls.args$setZE[which.max(abs(s))]
      return(list(c(1,rep(0,length(y)))))
    }
    
  }
  

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
  index.breakpoint <- which(duplicated(res.ode[,-(1:2)])) - 1# which(abs(res.ode[,2])>0) #
  
  if(length(index.breakpoint) == 0){ ## no breakpoint
    res <- list(breakpoint.sign = FALSE,
                breakpoint.constrain = FALSE,
                index_ode = n.lambda,
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
                 index_coef = c(setNE,setPE)[which.min(index.changeSign)],
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
      # H <- hessian(y) #### problem when using LAVA derivatives
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
                  index_coef = setZE[index],
                  sign_coef =  sign(s[index]),
                  newLambda = res.ode[index.breakpoint,1])
    }
  }
  
  #### post treatment
  res$index_plus <- min(n.lambda, res$index_ode + 1)
  res$index_minus <- max(1, res$index_ode - 1)
  
  #### export
  return(res)
  
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
