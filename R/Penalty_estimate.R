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
#' @param lambdaN Nuclear norm penalization parameter
#' @param regularizationPath the algorithm used to compute the regularization path. 
#' 0 indicates to estimate the solution for a given lambda
#' 1 corresponds to the algorithm proposed by (Park 2007) - called GLMpath.  Only works for regression models.  
#' 2 corresponds to the algorithm proposed by (Zhou 2014) - called EPSODE. 
#' If regularizationPath>0, the argument lambda1 is ignored but not lambda2
#' @param resolution_lambda1 argument for the EPSODE function (see Penalty_EPSODE.R)
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
#' EPSODE1lvm.fit <- estimate(plvm.model,  data = df.data, regularizationPath = 2)
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
estimate.plvm <- function(x, data, lambda1, lambda2, lambdaN, adaptive = FALSE, control = list(), estimator = "penalized", 
                          regularizationPath = FALSE, resolution_lambda1 = c(1e-1,1e-3), increasing = TRUE, stopLambda = NULL, stopParam = NULL, exportAllPath = FALSE, 
                          fit = "BIC",
                          method.proxGrad = "ISTA", step = 1, BT.n = 100, BT.eta = 0.8, force.descent = FALSE,
                          fixSigma = FALSE, ...) {
  
  names.coef <- coef(x)
  n.coef <- length(names.coef)
  
  #### prepare control
  if("iter.max" %in% names(control) == FALSE){
    control$iter.max <- 1000
  }
  if("trace" %in% names(control) == FALSE){
    control$trace <- 0
  }else{
    control$trace <- control$trace - 1
  }
  if("constrain" %in% names(control) == FALSE){
    control$constrain <- FALSE
  }
  
  #### prepare data (scaling)
  if(control$trace>=0){cat("Scale and center dataset \n")}
  resData <- prepareData.lvm(x, data = data, penalty = x$penalty)
  penalty <- resData$penalty
  
  ## penalty
  if(!missing(lambda1)){
    penalty$lambda1 <- as.numeric(lambda1)
  }
  if(!missing(lambda2)){
    penalty$lambda2 <- as.numeric(lambda2)
  }
  if(!missing(adaptive)){
    penalty$adaptive <- as.numeric(adaptive)
  }
  penalty$names.varCoef <- names.coef[x$index$parBelongsTo$cov]
  penalty$names.varCoef <- intersect(penalty$names.varCoef,
                                     paste(c(endogenous(x), latent(x)), 
                                           c(endogenous(x), latent(x)), 
                                           sep = ","))
  control$penalty <- penalty
  
  ## penaltyNuclear
  if(!is.null(x$penaltyNuclear$name.coef)){
    if(!missing(lambdaN)){
      x$penaltyNuclear$lambdaN <- as.numeric(lambdaN)
    }
    if(length(endogenous(x))>1){
      stop("nuclear norm only implemented for a linear model \n")
    }
    
    x$penaltyNuclear$objective <- function(coef){
      x$penaltyNuclear$FCTobjective(coef, 
                                    Y = as.vector(resData$data[,x$penaltyNuclear$name.Y,drop=TRUE]), 
                                    X = as.matrix(resData$data[,c(exogenous(x),x$penaltyNuclear$name.X),drop=FALSE]))}
    x$penaltyNuclear$gradient <- function(coef){
      x$penaltyNuclear$FCTgradient(coef, 
                                   Y = as.vector(resData$data[,x$penaltyNuclear$name.Y,drop=TRUE]), 
                                   X = as.matrix(resData$data[,c(exogenous(x),x$penaltyNuclear$name.X),drop=FALSE]))}
    control$penaltyNuclear <- x$penaltyNuclear
  }
  
  ## pass parameters for the penalty through the control argument
  control$proxGrad <- list(method = method.proxGrad,
                           step = step,
                           BT.n = BT.n,
                           BT.eta = BT.eta,
                           expX = if(control$constrain){penalty$names.varCoef}else{NULL},
                           force.descent = force.descent,
                           trace = if(regularizationPath == 0){control$trace}else{FALSE},
                           fixSigma = fixSigma,
                           envir = environment() # pass x and data in case where fixSigma = TRUE
  )
  
  control$regPath <- list(type = regularizationPath,
                          resolution_lambda1 = resolution_lambda1,
                          increasing = increasing,
                          stopParam = stopParam,
                          stopLambda = stopLambda,
                          fixSigma = fixSigma,
                          exportAllPath = exportAllPath,
                          trace = if(regularizationPath > 0){control$trace}else{FALSE})
  
  #### initialization
  if(is.null(control$start) || regularizationPath>0){
    if(control$trace>=0){cat("Initialization: ")}
    res.init  <- initializer.lvm(x, data = resData$data, names.coef = names.coef, n.coef = n.coef,
                                 penalty = control$penalty, regPath = control$regPath, ...)
    
    if(regularizationPath == 0){
      
      if(control$trace>=0){cat(" LVM where all penalized coefficient are shrinked to 0 \n")}
      control$start <- res.init$lambdaMax
      # if(control$trace>=0){cat(" unpenalized LVM \n")}
      # control$start <- res.init$lambdaMax
      
    }else {
      control$regPath$beta_lambdaMax <- res.init$lambdaMax
      control$regPath$beta_lambda0 <- res.init$lambda0
    }
  }
  
  ## proximal operator 
  resOperator <- init.proxOperator(lambda1 = control$penalty$lambda1,  # will be NULL if control$penalty does not exist
                                   lambda2 = control$penalty$lambda2, 
                                   lambdaN = control$penaltyNuclear$lambdaN,  # will be NULL if control$penaltyNuclear does not exist
                                   group.coef = control$penalty$group.coef,
                                   regularizationPath = regularizationPath)
  control$proxOperator <- resOperator$proxOperator
  control$objectivePenalty <- resOperator$objectivePenalty
  
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
    }
    
    return(res)
    
  }else{
    res <- lava:::estimate.lvm(x = x, data = resData$data, 
                               method = if(regularizationPath == 0){"optim.regLL"}else{"optim.regPath"}, 
                               control = control, estimator = estimator, ...)
  }
  class(res) <- append("plvmfit", class(res))
  
  ## add elements to object
  res$penalty <-  control$penalty[c("name.coef", "group.coef", "lambda1", "lambda2")]
  res$penalty$regularizationPath <- regularizationPath
  res$penalty$penaltyNuclear <- control$penaltyNuclear
  if( regularizationPath > 0){
    res$regularizationPath <- res$opt$message
    res$opt$message <- ""
  }else{
    if(fixSigma){
      res$penalty$lambda1.abs <- res$penalty$lambda1
      res$penalty$lambda2.abs <- res$penalty$lambda2
      res$penalty$lambda1 <- res$penalty$lambda1/sum(coef(res)[control$penalty$names.varCoef])
      res$penalty$lambda2 <- res$penalty$lambda2/sum(coef(res)[control$penalty$names.varCoef])
    }else{
      res$penalty$lambda1.abs <- res$penalty$lambda1*sum(coef(res)[control$penalty$names.varCoef])
      res$penalty$lambda2.abs <- res$penalty$lambda2*sum(coef(res)[control$penalty$names.varCoef])
    }
  }
  
  #### estimate the best model according to the fit parameter
  if(regularizationPath > 0 && !is.null(fit)){
    if(control$trace>=0){cat("Best penalized model according to the",fit,"criteria",if(control$trace>=1){"\n"})}
    res <- calcLambda(res, model = x, fit = fit, data.fit = data, trace = control$trace)
    if(control$trace>=0){cat(" - done \n")}
  }
  
  #### export
  res$x <- x
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
  if(nrow(data) > length(coef(x))){
    suppressWarnings(
      initLVM <- try(lava:::estimate.lvm(x = x, data = data, ...), silent = TRUE)
    )
    
    if(("try-error" %in% class(initLVM) == FALSE)){ # should also check convergence
      start_lambda0 <- coef(initLVM)
    }else{ 
      start_lambda0 <- NULL
    }
  }else{
    start_lambda0 <- NULL
  }
  
  #### hight dimensional model
  x0 <- x
  start_lambdaMax <- setNames(rep(0,n.coef),names.coef)
  
  ## removed penalized variables
  for(iter_link in penalty$name.coef){
    x0 <- rmLink(x0, iter_link)
  }
  
  ## estimate the model
  suppressWarnings(
    x0.fit <- lava:::estimate.lvm(x = x0, data = data, ...)
  )
  newCoef <- coef(x0.fit, level = 9)[,"Estimate"]
  newCoef <- newCoef[names(newCoef) %in% names(start_lambdaMax)] # only keep relevant parameters
  start_lambdaMax[names(newCoef)] <- newCoef
 # start_lambdaMax[] <- coef(x0.fit, level = 9)[,"Estimate"]
  
  return(list(lambda0 = start_lambda0,
              lambdaMax = start_lambdaMax))
}

#' @title Prepare the data for the estimate function
#'  
#' Scale the data, update the penalty term according to the presence of factors. If regularizationPath=1 (i.e. glmPath) orthogonolize the data
#' 
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param penalty parameter of the penalty
#' 
prepareData.lvm <- function(x, data, method.center = "mean", method.scale = "sd", penalty){
  
  #### rescale data
  varToScale <- names(which(unlist(lapply(data,is.numeric))))
  if(class(data)[1] != "data.frame"){data <- as.data.frame(data)}
  
  value.center <- sapply(varToScale, function(x){do.call(method.center,args = list(data[[x]]))})
  value.scale <- sapply(varToScale, function(x){do.call(method.scale,args = list(data[[x]]))})
  data[, varToScale] <- scale(data[, varToScale, drop = FALSE], center = value.center, scale = value.scale)
  # data[, c(varToScale) := lapply(.SD,scale), .SDcols = varToScale]
  
  #### handle factors for the penalty BUT DOES NOT WORK IF SEVERAL APPARTION OF THE SAME VARIABLE
  test.factor <- unlist(lapply(manifest(x), function(name){is.factor(data[[name]])}))
  if(any(test.factor)){ # if some manifest variables are factors
    for(iterFactor in manifest(x)[test.factor]){
      levelTempo <- levels(data[[iterFactor]])
      indexTempo <- grep(paste0("+~",iterFactor,"$"), penalty$name.coef)
      varTempo <- penalty$name.coef[indexTempo]
      
      if(length(varTempo)>0){ # if the manifest variable is penalized
        newTempo <- paste0(varTempo, levelTempo[-1]) # replace by the non reference level
        n.newTempo <- length(newTempo)
        
        penalty$V <- rbind(penalty$V,
                           matrix(0,nrow = n.newTempo, ncol = ncol(penalty$V))
        )
        penalty$V <- cbind(penalty$V,
                           matrix(0,nrow = nrow(penalty$V), ncol = n.newTempo)
        ) 
        colnames(penalty$V)[seq(ncol(penalty$V)-n.newTempo+1,ncol(penalty$V))] <- newTempo
        rownames(penalty$V)[seq(ncol(penalty$V)-n.newTempo+1,ncol(penalty$V))] <- newTempo
        penalty$V[,newTempo] <- penalty$V[,varTempo]
        penalty$V[newTempo,] <- penalty$V[,varTempo]
        for(iterNew in newTempo){
          penalty$V[iterNew,iterNew] <- penalty$V[varTempo,varTempo]
        }
        penalty$V <- penalty$V[rownames(penalty$V)!=varTempo,colnames(penalty$V)!=varTempo]
        
        penalty$name.coef <- c(setdiff(penalty$name.coef,
                                       varTempo),
                               newTempo) 
        
        penalty$group.coef <- c(penalty$group.coef[-indexTempo],
                                rep(penalty$group.coef[indexTempo], length(newTempo))
        )
        
      }
    }
  }
  
  #### Should use ???
  # if (length(exogenous(x) > 0)) {
  #   catx <- lava:::categorical2dummy(x, data)
  #   x <- catx$x
  #   data <- catx$data
  # }
  
  #### export
  return(list(data = data, 
              scale = value.scale,
              center = value.center,
              penalty = penalty))
}



