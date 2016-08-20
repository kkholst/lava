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
  penaltyNuclear <- control$penaltyNuclear
  control <- control[names(control) %in% PGcontrols]
  
  n.coef <- length(start)
  
  #### specify constrain: divide all variance parameters by the variance fo the first parameter and fix first variance parameter to one
  if(control$proxGrad$fixSigma){
    index.constrain <- which(names(start) %in% penalty$names.varCoef)
    if(control$trace>=0){cat("constrains: lambda1=lambda1/sum(", paste(names(start)[index.constrain], collapse = " "),")\n")}
  }else{
    index.constrain <- NULL
  }
  
  #### Proximal gradient algorithm
  # update penalty
  newPenalty <- initPenalty(start = start, penalty = penalty, penaltyNuclear = penaltyNuclear)
  
  proxOperator <- function(x, step){
    control$proxOperator(x, step = step,
                         lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1, lambda2 = newPenalty$lambda2, 
                         test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                         nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol,
                         index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
  }
  if(!is.na(newPenalty$lambdaN)){
    gradient <- penaltyNuclear$gradient
    objective <- penaltyNuclear$objective
    if(length(index.constrain)>0){index.constrain <- which(names(newPenalty$start) %in% names(start)[index.constrain])}
  }
  
  if(control$trace>=0){cat("Proximal gradient ")}
  res <- proxGrad(start = newPenalty$start, proxOperator = proxOperator, 
                  hessian = hessian, gradient = gradient, objective = objective,
                  step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta,
                  iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, force.descent = control$proxGrad$force.descent,
                  method = control$proxGrad$method, trace = control$proxGrad$trace)
  if(control$trace>=0){cat("- done \n")}
  
  if(penalty$adaptive){
    
    proxOperator <- function(x, step){
      control$proxOperator(x, step = step,
                           lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1/abs(res$par), lambda2 = newPenalty$lambda2, 
                           test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                           index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
    }
    
    if(control$trace>=0){cat("Proximal gradient (adaptive) ")}
    res <- proxGrad(start = res$par, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                    step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta,
                    iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, force.descent = control$proxGrad$force.descent,
                    method = control$proxGrad$method, trace = control$proxGrad$trace)
    if(control$trace>=0){cat("- done \n")}
    
  }
 
  ## estimation of the nuisance parameter
  if(length(index.constrain)>0){
    if(control$trace>=0){cat("Estimation of the nuisance parameter ")}
    res$par[index.constrain] <- optim.Nuisance.lvm(x = control$proxGrad$envir$x, data = control$proxGrad$envir$data, 
                                                   coefEstimated = res$par, coefNuisance = names(res$par)[index.constrain], coefPenalty = penalty$name.coef,
                                                   control = control, log = control$constrain)
    if(control$trace>=0){cat("- done \n")}
  }
  
  ## update objective with penalty - one value for each penalty
  objective.pen <- lapply(control$objectivePenalty, 
                          function(fct){do.call(fct, args = list(res$par, 
                                                                 lambdaN = penaltyNuclear$lambdaN, lambda1 = penalty$lambda1, lambda2 = penalty$lambda2, 
                                                                 test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2, test.penaltyN = newPenalty$test.penaltyN,
                                                                 nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol,
                                                                 expX = control$proxGrad$expX)
                          )})
  
  res$objective <- objective(res$par) + sum(unlist(objective.pen))
  
  ### export
  attr(res$message,"par") <- res$par
  res$par <- res$par[names(start)]
  return(res)
}

#' @title Estimate the penalization path
#' 
optim.regPath <- function(start, objective, gradient, hessian, control, ...){

  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  newPenalty <- initPenalty(start = start, penalty = control$penalty, regPath = control$regPath, penaltyNuclear = NULL)
 
    ##
  penalty <- newPenalty$penalty
  regPath <- newPenalty$regPath
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain", "proxGrad", "proxOperator", "regPath")
  control <- control[names(control) %in% PGcontrols]
 
  ## nuisance parameter
  if(regPath$fixSigma){
    if(length(penalty$names.varCoef)>1){stop("Cannot fixSigma when having several variance parameters \n")}
    indexNuisance <- which(names(start) %in% penalty$names.varCoef)
  }else{
    indexNuisance <- NULL
  }
 
  ## EPSODE
  resLassoPath <- EPSODE(beta_lambda0 = regPath$beta_lambda0, beta_lambdaMax = regPath$beta_lambdaMax,
                         objective = objective, gradient = gradient, hessian = hessian, 
                         V = penalty$V, 
                         indexPenalty = which(names(start) %in% penalty$name.coef), indexNuisance = indexNuisance, 
                         resolution_lambda1 = regPath$resolution_lambda1, increasing = regPath$increasing, stopLambda = regPath$stopLambda, stopParam = regPath$stopParam,
                         lambda2 = penalty$lambda2, test.penalty1 = newPenalty$test.penalty1,
                         control = control, exportAllPath = control$regPath$exportAllPath, trace = control$regPath$trace)

  
 ## estimation of the nuisance parameter and update lambda/lambda.abs
 if(length(indexNuisance)>0){
   if(control$trace>=0){cat("Estimation of the nuisance parameter ")}
   resLassoPath[,names(start)[indexNuisance]] <- sapply(1:NROW(resLassoPath), function(x){
     optim.Nuisance.lvm(x = control$proxGrad$envir$x, data = control$proxGrad$envir$data,
                        coefEstimated = unlist(resLassoPath[x,names(start)]),
                        coefNuisance = names(start)[indexNuisance], coefPenalty = penalty$name.coef,
                        control = control, log = FALSE)}
   )
   if(control$trace>=0){cat("- done \n")}
   resLassoPath$lambda1 <- resLassoPath$lambda1.abs/resLassoPath[,names(start)[indexNuisance]]
   resLassoPath$lambda2 <- resLassoPath$lambda2.abs/resLassoPath[,names(start)[indexNuisance]]
 }else{
   sumSigma <- apply(resLassoPath[,penalty$names.varCoef,drop=FALSE],1,sum)
   resLassoPath$lambda1.abs <- resLassoPath$lambda1*sumSigma
   resLassoPath$lambda2.abs <- resLassoPath$lambda2*sumSigma
 }
  
 #### export
  res <- list(par = start,
              convergence = 0,
              iterations = 0,
              evaluations = c("function" = 0, "gradient" = 0),
              message = resLassoPath[order(resLassoPath$lambda1),,drop = FALSE],
              objective = NA)
  
  ### export
  return(res)
}

#' @title Estimate nuisance parameter
optim.Nuisance.lvm <- function(x, data, 
                              coefEstimated, coefNuisance, coefPenalty,
                              control, log){
  
  control$trace <- FALSE
  names.coef <- setdiff(names(coefEstimated), coefNuisance)
  n.coef <- length(names.coef)
  
  fit.coef <- coefEstimated[names.coef, drop = FALSE]
 
  
  xConstrain <- x
  class(xConstrain) <- "lvm"
  
  for(iter_p in 1:n.coef){
    # if(fit.coef[iter_p]==0){
    #   xConstrain <- rmLink(xConstrain, var1 = names.coef[iter_p])
    # }else{
      xConstrain <- setLink(xConstrain, var1 = names.coef[iter_p], value = fit.coef[iter_p])
    # }
  }
  elvm <- estimate(xConstrain, data = data, control = control)
  
  if(log){
    return(log(coef(elvm)[coefNuisance]))
  }else{
    return(coef(elvm)[coefNuisance])
  }
  
}

#' @title initialise the penalty 
#' 
initPenalty <- function(start, regPath = NULL, penalty, penaltyNuclear){
  
  
  ## lasso
  if(!is.null(penalty$name.coef)){
    n.coef <- length(start)
    
    ## check
    if(any(penalty$name.coef %in% names(start) == FALSE)){
      warning("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
              "non-applied penalty: ",paste(setdiff(penalty$name.coef, names(start)), collapse = " "),"\n")
      penalty$group.coef <- penalty$group.coef[penalty$name.coef %in% names(start)]
      penalty$name.coef <- penalty$name.coef[penalty$name.coef %in% names(start)]
     }
   
    if(!is.null(regPath)){
      penalty$V <- penalty$V[names(start), names(start),drop = FALSE]
      regPath$beta_lambda0 <- regPath$beta_lambda0[names(start), drop = FALSE]
      regPath$beta_lambdaMax <- regPath$beta_lambdaMax[names(start), drop = FALSE]
    }
    
    lambdaN <- NA
    lambda1 <- rep(0, n.coef)
    lambda1[names(start) %in% penalty$name.coef] <- penalty$lambda1
    lambda2 <- rep(0, n.coef)
    lambda2[names(start) %in% penalty$name.coef] <- penalty$lambda2
    
    ## grouped lasso: set lasso indexes to 0
    test.penaltyN <- NULL
    test.penalty1 <- setNames(rep(0, n.coef), names(start))
    test.penalty1[penalty$name.coef] <- penalty$group.coef
    test.penalty2 <- lambda2>0
  }
  
  if(!is.null(penaltyNuclear$name.coef)){
    
    if(any(penaltyNuclear$name.coef %in% names(start) == FALSE)){
    start <- c(start[setdiff(names(start),penalty$names.varCoef)],
               setNames(rep(0,length(penaltyNuclear$name.coef)),penaltyNuclear$name.coef),
               start[penalty$names.varCoef])  
    }
    
    n.coef <- length(start)
    
    lambdaN <- penaltyNuclear$lambdaN
    lambda1 <- NA
    lambda2 <- NA
   
    test.penaltyN <- which(names(start) %in% penaltyNuclear$name.coef)
    test.penalty1 <- NULL
    test.penalty2 <- NULL
  }
  
  return(list(start = start, regPath = regPath, penalty = penalty,
              lambdaN = lambdaN, lambda1 = lambda1, lambda2 = lambda2,
              test.penaltyN = test.penaltyN, test.penalty1 = test.penalty1, test.penalty2 = test.penalty2))
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
