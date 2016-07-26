#### OVERVIEW
# optim.proxGrad: Estimate a penalized lvm model using a proximal gradient algorithm 
# optim.regPath: Estimate the regularization path associated to a LVM
# optim.nuisance: Estimate a nuisance parameter given all others
# initPenalty: initialise the penalty 

#' @title Estimate a penalized lvm model using a proximal gradient algorithm
#' 
optim.regLL <- function(start, objective, gradient, hessian, control, ...){
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain", "proxOperator", "objectivePenalty",
                  "proxGrad")
  
  penalty <- control$penalty
  control <- control[names(control) %in% PGcontrols]
  
  n.coef <- length(start)
  
  #### specify constrain: divide all variance parameters by the variance fo the first parameter and fix first variance parameter to one
  if(control$proxGrad$fixSigma){
    constrain <- which(names(start) %in% penalty$names.varCoef)
    if(control$trace>=0){cat("constrains: lambda1=lambda1/sum(", paste(names(start)[constrain], collapse = " "),")\n")}
  }else{
    constrain <- NULL
  }
  
  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  res <- initPenalty(start = start, penalty = penalty)
  penalty$group.penaltyCoef <- res$group.penaltyCoef
  index.penaltyCoef <- res$index.penaltyCoef
  
  #### Proximal gradient algorithm
  lambda1 <- rep(0, n.coef)
  lambda1[index.penaltyCoef] <- penalty$lambda1
  lambda2 <- rep(0, n.coef)
  lambda2[index.penaltyCoef] <- penalty$lambda2
  
  test.penalty1<- penalty$group.penaltyCoef
  test.penalty2<- lambda2>0
  
 
  
  proxOperator <- function(x, step){ 
    if(length(constrain)>0){
      if(control$constrain){norm <- sum(exp(x[constrain]))
      }else{
        norm <- sum(x[constrain])
        if(norm < 0){stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")}
      }
    }else{
      norm <- 1
    }
    control$proxOperator(x, step,
                         lambda1 = lambda1/norm, lambda2 = lambda2/norm, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2, expX = control$proxGrad$expX)
  }
  
  res <- proxGrad(start = start, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                  step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, 
                  iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, force.descent = control$proxGrad$force.descent, 
                  method = control$proxGrad$method, trace = control$proxGrad$trace)
  
  
  if(penalty$adaptive){
    proxOperator <- function(x, step){
      if(length(constrain)>0){
        if(control$constrain){norm <- sum(exp(x[constrain]))
        }else{
          norm <- sum(x[constrain])
          if(norm < 0){stop("proxGrad: negative variance parameter - set constrain to TRUE in control \n")}
        }
      }else{
        norm <- 1
      }
      control$proxOperator(x, step,
                           lambda1 = lambda1/(norm*abs(res$par)), lambda2 = lambda2PG/norm, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2, expX = control$proxGrad$expX)
    }

    res <- proxGrad(start = res$par, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                    step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta,
                    iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, force.descent = control$proxGrad$force.descent,
                    method = control$proxGrad$method, trace = control$proxGrad$trace)

  }
 
  res$objective <- objective(res$par) + control$objectivePenalty(res$par, lambda1 = lambda1, lambda2 = lambda2, 
                                                                 test.penalty1 = penalty$group.penaltyCoef, test.penalty2 = lambda2>0, expX = control$proxGrad$expX)
  ### export
  return(res)
}

#' @title Estimate the penalization path
#' 
optim.regPath <- function(start, objective, gradient, hessian, control, ...){
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain", "proxGrad", "proxOperator", "regPath")
  
  penalty <- control$penalty
  regPath <- control$regPath
  control <- control[names(control) %in% PGcontrols]
  
  n.coef <- length(start)
  regPath$beta_lambdaMax <- control$regPath$beta_lambdaMax[names(start)] # remove the parameters fixed for identifiability purpose
  
  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  res <- initPenalty(start = start, penalty = penalty)
  penalty$group.penaltyCoef <- res$group.penaltyCoef
  index.penaltyCoef <- res$index.penaltyCoef
  
  #### main
  if( regPath$type == 1){
    resLassoPath <- glmPath(beta0 = start, objective = objective, gradient = gradient, hessian = hessian,
                            indexPenalty = index.penaltyCoef, indexNuisance = which(names(start) %in% penalty$names.varCoef), 
                            sd.X = penalty$sd.X, base.lambda1 = penalty$scaleLambda1, lambda2 = penalty$scaleLambda2*penalty$lambda2, group.lambda1 = penalty$group.penaltyCoef,
                            control = control)
    
  }else if( regPath$type == 2){
    V <- penalty$V[names(start),names(start), drop = FALSE]
    group.penaltyCoef <- penalty$group.penaltyCoef[names(start)]
   
    ## check lambda
    if(regPath$increasing == TRUE){
      if(any(start[index.penaltyCoef] == 0)){
        stop("All penalized coefficient should be non-0 when using increasing = TRUE \n",
             "May be due to automatic initialization - in such a case set increasing to FALSE \n")
      }
    }else{
      if(!all(start[index.penaltyCoef] == 0)){
        stop("All penalized coefficient should be set to 0 when using increasing = FALSE \n")
      }
    }
    
    ## nuisance parameter
    if(is.null(regPath$fixSigma)){
      if(length(penalty$names.varCoef)==1){fixSigma <- TRUE}else{fixSigma <- FALSE}
    }
    
    if(regPath$fixSigma){
      indexNuisance <- which(names(start) %in% penalty$names.varCoef)
    }else{
      indexNuisance <- NULL
    }
    
    ## EPSODE
    resLassoPath <- EPSODE(beta = start, beta_lambdaMax = regPath$beta_lambdaMax, objective = objective, gradient = gradient, hessian = hessian, 
                           V = V, 
                           indexPenalty = index.penaltyCoef, indexNuisance = indexNuisance, 
                           stepLambda1 = regPath$stepLambda1, increasing = regPath$increasing, stopLambda = regPath$stopLambda, stopParam = regPath$stopParam,
                           lambda2 = penalty$lambda2, group.lambda1 = group.penaltyCoef,
                           control = control, trace = control$regPath$trace)
    
  }
 
  #### export
  res <- list(par = start,
              convergence = 0,
              iterations = 0,
              evaluations = c("function" = 0, "gradient" = 0),
              message = resLassoPath[order(resLassoPath$lambda1.abs),,drop = FALSE],
              objective = NA)
  
  ### export
  return(res)
}

#' @title initialise the penalty 
#' 
initPenalty <- function(start, penalty){
  
  n.coef <- length(start)
  
  ## check
  if(length(setdiff(penalty$names.penaltyCoef, names(start)))>0){
    warning("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
            "non-applied penalty: ",paste(setdiff(penalty$names.penaltyCoef, names(start)), collapse = " "),"\n")
    penalty$group.penaltyCoef <- penalty$group.penaltyCoef[penalty$names.penaltyCoef %in% names(start)]
    penalty$names.penaltyCoef <- penalty$names.penaltyCoef[penalty$names.penaltyCoef %in% names(start)]
  }
  index.penaltyCoef <- which(names(start) %in% penalty$names.penaltyCoef)
  
  ## grouped lasso: set lasso indexes to 0
  group.penaltyCoef <- setNames(rep(0, n.coef), names(start))
  group.penaltyCoef[penalty$names.penaltyCoef] <- penalty$group.penaltyCoef
  
  return(list(index.penaltyCoef = index.penaltyCoef, group.penaltyCoef = group.penaltyCoef))
}

initSigmaConstrain <- function(start, constrain, indexNuisance){
  
  if(constrain){
    start[indexNuisance] <- start[indexNuisance] - start[indexNuisance[1]]
    constrain <- setNames(0, names(start)[indexNuisance[1]]) 
  }else{
    start[indexNuisance] <- start[indexNuisance]/start[indexNuisance[1]]
    constrain <- setNames(1, names(start)[indexNuisance[1]]) 
  }
  indexAllCoef <- setdiff(1:length(start), indexNuisance[1])
  
  return(list(start = start,
              constrain = constrain,
              indexAllCoef = indexAllCoef)
  )
}
