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
  
  res <- proxGrad(start = newPenalty$start, proxOperator = proxOperator, 
                  hessian = hessian, gradient = gradient, objective = objective,
                  step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta,
                  iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, force.descent = control$proxGrad$force.descent,
                  method = control$proxGrad$method, trace = control$proxGrad$trace)

  
  # ### TO BE REMOVED
  # start2 <- c(start, setNames(rep(0, length(penaltyNuclear$name.coef)), penaltyNuclear$name.coef))
  # test.penaltyN <- names(start2) %in% penaltyNuclear$name.X
  # 
  # proxOperator <-  function(x, step){
  #   x[test.penaltyN] <- proxNuclear(x = x[test.penaltyN], step = step,
  #                                   lambda = penaltyNuclear$lambdaN, 
  #                                   nrow = penaltyNuclear$nrow, ncol = penaltyNuclear$ncol)
  #   return(x)
  # }
  # 
  # penaltyNuclear$objective(start2)
  # 
  # browser()
  # resLV <- proxGrad(start = start2, proxOperator,
  #                   objective = penaltyNuclear$objective, gradient = penaltyNuclear$gradient, 
  #                   step = 1, BT.n = 200, BT.eta = 0.5, force.descent = TRUE,
  #                   iter.max = 10, abs.tol = 1e-9, rel.tol = 1e-10, method = "ISTA", trace = TRUE)
  # 
  # 
  # ###
  
  
  
  if(penalty$adaptive){
    
    proxOperator <- function(x, step){
      control$proxOperator(x, step = step,
                           lambdaN = newPenalty$lambdaN, lambda1 = newPenalty$lambda1/abs(res$par), lambda2 = newPenalty$lambda2, 
                           test.penaltyN = newPenalty$test.penaltyN, test.penalty1 = newPenalty$test.penalty1, test.penalty2 = newPenalty$test.penalty2,
                           index.constrain = index.constrain, type.constrain = control$constrain, expX = control$proxGrad$expX)
    }
    
    res <- proxGrad(start = res$par, proxOperator = proxOperator, hessian = hessian, gradient = gradient, objective = objective,
                    step = control$proxGrad$step, BT.n = control$proxGrad$BT.n, BT.eta = control$proxGrad$BT.eta,
                    iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol, force.descent = control$proxGrad$force.descent,
                    method = control$proxGrad$method, trace = control$proxGrad$trace)

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
  
  PGcontrols <- c("iter.max","trace","abs.tol","rel.tol", "constrain", "proxGrad", "proxOperator", "regPath")
  
  penalty <- control$penalty
  regPath <- control$regPath
  control <- control[names(control) %in% PGcontrols]
  
  n.coef <- length(start)
  regPath$beta_lambdaMax <- control$regPath$beta_lambdaMax[names(start)] # remove the parameters fixed for identifiability purpose
  
  #### update the penalty according to start 
  # (some coefficient have been removed as they are chosen as a reference)
  newPenalty <- initPenalty(start = start, penalty = penalty, penaltyNuclear = NULL)
  
  #### main
  if( regPath$type == 1){
    resLassoPath <- glmPath(beta0 = start, objective = objective, gradient = gradient, hessian = hessian,
                            indexPenalty = which(names(start) %in% penalty$name.coef), indexNuisance = which(names(start) %in% penalty$names.varCoef), 
                            sd.X = penalty$sd.X, base.lambda1 = penalty$scaleLambda1, lambda2 = penalty$scaleLambda2*penalty$lambda2, test.penalty1 = newPenalty$test.penalty1,
                            control = control)
    
  }else if( regPath$type == 2){
    V <- penalty$V[names(start),names(start), drop = FALSE]
    group.coef <- penalty$group.coef[names(start)]
   
    ## check lambda
    if(regPath$increasing == TRUE){
      if(any(start[names(start) %in% penalty$name.coef] == 0)){
        stop("All penalized coefficient should be non-0 when using increasing = TRUE \n",
             "May be due to automatic initialization - in such a case set increasing to FALSE \n")
      }
    }else{
      if(!all(start[names(start) %in% penalty$name.coef] == 0)){
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
    resLassoPath <- EPSODE(beta = start, beta_lambdaMax = regPath$beta_lambdaMax,  beta_lambda0 = regPath$beta_lambda0,
                           objective = objective, gradient = gradient, hessian = hessian, 
                           V = V, 
                           indexPenalty = which(names(start) %in% penalty$name.coef), indexNuisance = indexNuisance, 
                           resolution_lambda1 = regPath$resolution_lambda1, increasing = regPath$increasing, stopLambda = regPath$stopLambda, stopParam = regPath$stopParam,
                           lambda2 = penalty$lambda2, test.penalty1 = newPenalty$test.penalty1,
                           control = control, exportAllPath = control$regPath$exportAllPath, trace = control$regPath$trace)
    
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
initPenalty <- function(start, penalty, penaltyNuclear){
  
  
  ## lasso
  if(!is.null(penalty$name.coef)){
    n.coef <- length(start)
    
    ## check
    if(length(setdiff(penalty$name.coef, names(start)))>0){
      warning("initPenalty: some penalty will not be applied because the corresponding parameter is used as a reference \n",
              "non-applied penalty: ",paste(setdiff(penalty$name.coef, names(start)), collapse = " "),"\n")
      penalty$group.coef <- penalty$group.coef[penalty$name.coef %in% names(start)]
      penalty$name.coef <- penalty$name.coef[penalty$name.coef %in% names(start)]
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
    start <- c(start[setdiff(names(start),penalty$names.varCoef)],
               setNames(rep(0,length(penaltyNuclear$name.coef)),penaltyNuclear$name.coef),
               start[penalty$names.varCoef])  
    
    n.coef <- length(start)
    
    lambdaN <- penaltyNuclear$lambdaN
    lambda1 <- NA
    lambda2 <- NA
   
    test.penaltyN <- which(names(start) %in% penaltyNuclear$name.coef)
    test.penalty1 <- NULL
    test.penalty2 <- NULL
  }
  
  return(list(start = start, 
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
