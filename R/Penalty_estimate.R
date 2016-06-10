#### OVERVIEW
# estimate.plvm: Estimate a penalized lvm model
# initializer.lvm: if not given by the user, find an initial solution to LVM that will be used to initialize the regularization path algorithm or the proximal gradient algorithm
# orthoData.lvm: [ONLY used if regularizationPath = 1] make the non-penalized variables orthogonal to the penalized variables
# rescaleRes: [ONLY used if regularizationPath = 1]  Cancel the effect of the orthogonalization on the estimated parameters

#' @title Estimate a penalized lvm model
#
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param lambda1 L1 penalization parameter
#' @param lambda2 L2 penalization parameter
#' @param regularizationPath the algorithm used to compute the regularization path. 
#' 0 indicates to estimate the solution for a given lambda
#' 1 corresponds to the algorithm proposed by (Park 2007) - called GLMpath.  Only works for regression models.  
#' 2 corresponds to the algorithm proposed by (Zhou 2014) - called EPSODE. 
#' If regularizationPath>0, the argument lambda1 is ignored but not lambda2
#' @param stepLambda1 argument for the EPSODE function (see Penalty_EPSODE.R)
#' @param method.proxGrad argument for the proxGrad function (see Penalty_optims.R) 
#' @param step argument for the proxGrad function (see Penalty_optims.R)
#' @param BT.n argument for the proxGrad function (see Penalty_optims.R)
#' @param fixSigma should the variance parameter be fixed at 1 ? Only works for regression models. [temporary]
#' @param ... additional arguments to be passed to lava:::estimate.lvm
#' 
#' @details 
#' 
#' @references 
#' Zhou 2014 - A generic Path Algorithm for Regularized Statistical Estimation
#' Park 2007 - L1-regularization path algorithm for GLM
#' 
#' @return 
#' a plvmfit object
#' 
#' @examples
#' set.seed(10)
#' n <- 300
#' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:4), collapse = "+")))
#' lvm.modelSim <- lvm(formula.lvm)
#' df.data <- sim(lvm.modelSim,n)
#' 
#' lvm.model <- lvm(formula.lvm)
#' plvm.model <- penalize(lvm.model)
#' 
#' #### unpenalized
#' lvm.fit <- estimate(lvm.model, df.data)
#' 
#' #### L1 penalisation
#' 
#' ## regularization path
#' rp1lvm.fit <- estimate(plvm.model,  data = df.data, regularizationPath = 1)
#' rp1lvm.fit
#' 
#' EPSODE1lvm.fit <- estimate(plvm.model,  data = df.data, stepLambda1 = 100, regularizationPath = 2)
#' EPSODE1lvm.fit
#' 
#' ## 
#' p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1)
#' p1lvm.fit
#' 
#' p1lvm.fit_bis <- estimate(plvm.model,  data = df.data, lambda1 = rp1lvm.fit$opt$message[2,"lambda1.abs"])
#' p1lvm.fit_bis
#' 
#' #### L2 penalisation
#' p2lvm.fit <- estimate(plvm.model,  data = df.data, lambda2 = 0.1)
#' p2lvm.fit
#' 
#' #### elastic net
#' rp12lvm.fit <- estimate(plvm.model,  data = df.data, regularizationPath = 1, lambda2 = 5, fixSigma = TRUE)
#' rp12lvm.fit
#' 
#' p12lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1, lambda2 = 0.1)
#' p12lvm.fit
estimate.plvm <- function(x, data, lambda1, lambda2, control = list(),
                          regularizationPath = FALSE, stepLambda1 = NULL, increasing = TRUE,
                          method.proxGrad = "ISTA", step = 1, BT.n = 20, BT.eta = 0.5, trace = FALSE,
                          fixSigma = FALSE, ...) {
   
  names.coef <- coef(x)
  n.coef <- length(names.coef)
  
  #### prepare arguments
  if("iter.max" %in% names(control) == FALSE){
    control$iter.max <- 1000
  }
  
  ## penalty
  penalty  <- x$penalty
  
  if(!missing(lambda1)){
    penalty$lambda1 <- as.numeric(lambda1)
  }
  if(!missing(lambda2)){
    penalty$lambda2 <- as.numeric(lambda2)
  }
  penalty$names.varCoef <- names.coef[x$index$parBelongsTo$cov]
  
  if(regularizationPath > 0){
    save_lambda2 <-  penalty$lambda2
  }
  
  control$penalty <- penalty
  
  ## pass parameters for the penalty through the control argument
  control$proxGrad <- list(method = method.proxGrad,
                           step = step,
                           BT.n = BT.n,
                           BT.eta = BT.eta,
                           trace = if(regularizationPath == 0){trace}else{FALSE},
                           fixSigma = fixSigma
  )
  
  control$regPath <- list(type = regularizationPath,
                          stepLambda1 = stepLambda1,
                          increasing = increasing,
                          trace = if(regularizationPath > 0){trace}else{FALSE})
  
  ## proximal operator 
  resOperator <- init.proxOperator(lambda1 = penalty$lambda1, 
                                   lambda2 = penalty$lambda2, 
                                   group.penaltyCoef = penalty$group.penaltyCoef,
                                   regularizationPath = regularizationPath)
  control$proxOperator <- resOperator$proxOperator
  control$objectivePenalty <- resOperator$objectivePenalty
  
  #### orthogonolize data for glmPath
  if(regularizationPath == 1){
    resOrtho <- prepareDataPath.lvm(model = x, data = data, penalty = penalty)
    
    data <- resOrtho$data
    control$start <- resOrtho$mu.X
    x <- resOrtho$model
    control$penalty <- resOrtho$penalty
    control$penalty$lambda2 <- control$penalty$lambda2 * save_lambda2
    names.coef <- coef(x)
  }
  
  #### initialization
  if(is.null(control$start)){
    res.init  <- initializer.lvm(x, data = data, names.coef = names.coef, n.coef = n.coef,
                                 penalty = control$penalty, regPath = control$regPath, ...)
    
    if(regularizationPath > 0){
      control$start <- res.init$lambda0
      
      if(regularizationPath == 2){
        control$regPath$beta_lambdaMax <- res.init$lambdaMax
      }
    }else{
      control$start <- res.init$lambdaMax
    }
   
  }

  #### main
  if(all(c("objective", "gradient", "hessian") %in% names(list(...)))){
    
    if(regularizationPath == 0){
      res <- optim.regLL(start = control$start, 
                            objective = list(...)$objective, gradient = list(...)$gradient, hessian = list(...)$hessian, 
                            control = control)
    }else{
      res <- optim.regPath(start = control$start, 
                           objective = list(...)$objective, gradient = list(...)$gradient, hessian = list(...)$hessian, 
                           control = control)
      
      if(regularizationPath == 1){
        res <- rescaleRes(Mres = as.matrix(res), 
                          penalty = control$penalty, 
                          orthogonalizer = resOrtho$orthogonalizer)
      }
    }
    
    return(res)
    
  }else{
    res <- lava:::estimate.lvm(x = x, data = data, estimator = "penalized", 
                               method = if(regularizationPath == 0){"optim.regLL"}else{"optim.regPath"}, 
                               control = control, ...)
  }
  class(res) <- append("plvmfit", class(res))
  
  #### update sigma value
  if(regularizationPath > 0 || fixSigma == TRUE){
    res <- estimate.nuisance(x, plvmfit = res, data = data, control = control, 
                             regularizationPath = regularizationPath)
    
  }
  
  #### rescale parameter to remove the effect of orthogonalization
  if(regularizationPath == 1){
      res$opt$message <- rescaleRes(Mres = as.matrix(res$opt$message), 
                                  penalty = control$penalty, 
                                  orthogonalizer = resOrtho$orthogonalizer)
  }
  
  #### export
  res$penalty <-  control$penalty[c("names.penaltyCoef", "group.penaltyCoef", "lambda1", "lambda2")]
  res$penalty$regularizationPath <- regularizationPath
 
  return(res)
}


#' @title Find an initial solution for the model
#
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param n.coef number of coefficients
#' @param penalty parameter of the penalty
#' @param regPath parameter of the penalisation path
#' 
initializer.lvm <- function(x, data, names.coef, n.coef, penalty, regPath, ...){

  #### normal model
  initLVM <- try(lava:::estimate.lvm(x = x, data = data, ...))
  if("try-error" %in% class(initLVM) == FALSE){
    start_lambda0 <- coef(initLVM)
  }else{ 
    start_lambda0 <- NULL
  }
    
  #### hight dimensional model
  x0 <- x
  start_lambdaMax <- setNames(rep(0,n.coef),names.coef)
  
  ## removed penalized variables
  for(iter_link in penalty$names.penaltyCoef){
    
    if(iter_link %in% coef(x)[x$index$parBelongsTo$cov] ){
      kill(x0) <- strsplit(iter_link, split = ",", fixed = TRUE)[[1]][1]
      start[iter_link] <- 1
    }else if(iter_link %in% coef(x)[x$index$parBelongsTo$mean]){
      cancel(x0) <- as.formula(paste0(iter_link,"~1"))
    }else {
      cancel(x0) <- as.formula(iter_link)
    }
  }
  
  ## estimate the model
  suppressWarnings(
    x0.fit <- lava:::estimate.lvm(x = x0, data = data, ...)
  )
  start_lambdaMax[names(coef(x0.fit))] <- coef(x0.fit)

  return(list(lambda0 = start_lambda0,
              lambdaMax = start_lambdaMax))
}


#' @title Estimate nuisance parameter
#' 
estimate.nuisance <- function(x, plvmfit, data, control, regularizationPath){
  
  xConstrain <- x
  class(xConstrain) <- "lvm"
  if(regularizationPath){
    fit.coef <- penPath(plvmfit, type = "coef")
    names.coef <- names(fit.coef)
  }else{
    names.coef <- names( coef(plvmfit))
    fit.coef <- setNames(data.frame(rbind(coef(plvmfit))), names.coef)
  }
  
  
  coefReg <- fit.coef[,grep("~", names.coef, fixed = TRUE), drop = FALSE]
  names.coefReg <- names(coefReg)
  coefInt <- fit.coef[,setdiff(1:length(names.coef), grep("~|,", names.coef)), drop = FALSE]
  names.coefInt <- names(coefInt)
  names.coefVar <- setdiff(names.coef, c(names.coefReg,names.coefInt))
  
  for(iterPath in 1:nrow(fit.coef)){
    
    for(iter_p in 1:length(coefReg)){
      
      regression(xConstrain, as.formula(names.coefReg[iter_p])) <- coefReg[iterPath,iter_p]
    }
    
    for(iter_p in 1:length(coefInt)){
      intercept(xConstrain, as.formula(paste0("~",names.coefInt[iter_p]))) <- coefInt[iterPath, iter_p]
    }
    
    plvmfit2 <- estimate(xConstrain, data = data, control = control)
    
    if(regularizationPath){
      penPath(plvmfit, row = iterPath) <- coef(plvmfit2)[names.coefVar]
      penPath(plvmfit, names = "lambda1", row = iterPath) <- penPath(plvmfit, row = iterPath, type = "lambda1.abs") / fit.coef[iterPath,names.coefVar[1]]
      penPath(plvmfit, names = "lambda2", row = iterPath) <- penPath(plvmfit, row = iterPath, type = "lambda2.abs") / fit.coef[iterPath,names.coefVar[1]]
    }else{
      index1 <- match(names(coef(plvmfit)), names.coefVar, nomatch = 0)>0
      index2 <- match(names(coef(plvmfit2)), names.coefVar, nomatch = 0)>0
      plvmfit$coef[index1,] <- plvmfit2$coef[index2,]
    }
  }
  
  return(plvmfit)
  
}
