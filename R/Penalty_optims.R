#### OVERVIEW
# optim.proxGrad: Estimate a penalized lvm model using a proximal gradient algorithm 
# optim.regPath: Estimate the regularization path associated to a LVM
# optim.nuisance: Estimate a nuisance parameter given all others
# initPenalty: initialise the penalty 

#' @title Estimate a penalized lvm model using a proximal gradient algorithm
#' 
optim.proxGrad <- function(start, objective, gradient, hessian,...){
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain", "proxOperator", "proxGrad")
 
  dots <- list(...)
  penalty <- dots$control$penalty
  control <- dots$control[names(dots$control) %in% PGcontrols]
  
  n.coef <- length(start)
  
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
  
  res <- ISTA(start = start, proxOperator = control$proxOperator, hessian = hessian, gradient = gradient, objective = objective,
              lambda1 = lambda1, lambda2 = lambda2, group.lambda1 = penalty$group.penaltyCoef, constrain = control$proxGrad$fixSigma,
              step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta, trace = control$proxGrad$trace,
              iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, fast = control$proxGrad$fast)
 
  res$objective <- objective(res$par) + penalized_objectivePen.lvm(res$par, lambda1 = lambda1, lambda2 = lambda2)
  
    ## estimate the variance parameter if it has been fixed before
  if(!is.null(control$proxGrad$fixSigma)){
    index.sigma2 <- which(names(res$par) == names(control$proxGrad$fixSigma)) 
 
    res2 <- optim.nuisance(coef = res$par, 
                          index.sigma2 = index.sigma2,
                          sigma_max = if(control$constrain){NULL}else{control$proxGrad$sigmaMax},
                          objective = objective,
                          gradient = gradient)

    res$par[index.sigma2] <- res2
  }
  
  ### export
  return(res)
}

#' @title Estimate the penalization path
#' 
optim.regPath <- function(start, objective, gradient, hessian,...){
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain", "proxGrad", "proxOperator", "regPath")
  
  dots <- list(...)
  penalty <- dots$control$penalty
  regPath <- dots$control$regPath
  control <- dots$control[names(dots$control) %in% PGcontrols]
  
  n.coef <- length(start)

  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  res <- initPenalty(start = start, penalty = penalty)
  penalty$group.penaltyCoef <- res$group.penaltyCoef
  
  index.penaltyCoef <- res$index.penaltyCoef
  indexNuisance <- which(names(start) == names(control$proxGrad$fixSigma))
  
  #### main
  if( regPath$type == 1){
    resLassoPath <- glmPath(beta0 = start, objective = objective, gradient = gradient, hessian = hessian,
                            indexPenalty = index.penaltyCoef, indexNuisance = indexNuisance, 
                            sd.X = penalty$sd.X, base.lambda1 = penalty$lambda1, lambda2 = penalty$lambda2, group.lambda1 = penalty$group.penaltyCoef,
                            control = control)
    
  }else if( regPath$type == 2){
    
    V <- matrix(0, nrow = length(start), ncol = length(start))
    diag(V)[index.penaltyCoef] <- 1
    
    resLassoPath <- EPSODE(beta = start, objective = objective, gradient = gradient, hessian = hessian, V = V, 
                           indexPenalty = index.penaltyCoef, indexNuisance = indexNuisance, 
                           stepLambda1 = regPath$stepLambda1, correctionStep = regPath$correctionStep,
                           lambda2 = penalty$lambda2, group.lambda1 = penalty$group.penaltyCoef,
                           control = control)
    
  }
  
  #### update sigma value
  if(!is.null(control$proxGrad$fixSigma)){
  resLassoPath[,names(control$proxGrad$fixSigma)] <- apply(resLassoPath[,names(start), drop = FALSE], 1, optim.nuisance, 
                                                           index.sigma2 = indexNuisance,
                                                           sigma_max = if(control$constrain){NULL}else{control$proxGrad$sigmaMax},
                                                           gradient = gradient,
                                                           objective = objective)
  
  if(control$constrain){
    resLassoPath[,names(control$proxGrad$fixSigma)] <- exp(resLassoPath[,names(control$proxGrad$fixSigma)])
  }
  
  resLassoPath$lambda1 <- resLassoPath$lambda1.abs / resLassoPath[,names(control$proxGrad$fixSigma)]
  }
  
 
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
#' 
optim.nuisance <- function(coef, index.sigma2, sigma_max, objective, gradient, method = "L-BFGS-B"){

  warperObj <- function(sigma2){
    pp <- coef
    pp[index.sigma2] <- sigma2
    return(objective(pp))
  }
  
  warperGrad <- function(sigma2){
    pp <- coef
    pp[index.sigma2] <- sigma2
    return(gradient(pp)[index.sigma2])
  }
  
  if(is.null(sigma_max)){
    coef[index.sigma2] <- optim(par = coef[index.sigma2], 
                                fn = warperObj, gr = warperGrad, method = method)$par
  }else{
    coef[index.sigma2] <- optim(par = if(coef[index.sigma2]>=sigma_max){sigma_max/2}else{coef[index.sigma2]}, 
                                fn = warperObj, gr = warperGrad, lower = 1e-12, upper = sigma_max, method = method)$par
  }
    

#   coef[index.sigma2] <- optim(par = coef[index.sigma2], 
#                               fn = warperGrad, method = method,
#                               lower = 1e-12)$par # , upper = sigma_max

}

#' @title initialise the penalty 
#' 
initPenalty <- function(start, penalty){
  
  n.coef <- length(start)
  
  ## check
  if(length(setdiff(penalty$names.penaltyCoef, names(start)))>0){
    message("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
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


