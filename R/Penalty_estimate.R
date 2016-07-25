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
estimate.plvm <- function(x, data, lambda1, lambda2, adaptive = FALSE, control = list(),
                          regularizationPath = FALSE, stepLambda1 = NULL, increasing = TRUE, stopLambda = NULL, stopParam = NULL, fit = "BIC",
                          method.proxGrad = "ISTA", step = 1, BT.n = 100, BT.eta = 0.8, force.descent = FALSE,
                          fixSigma = FALSE, ...) {
  
  names.coef <- coef(x)
  n.coef <- length(names.coef)
  
  #### prepare control
  if("iter.max" %in% names(control) == FALSE){
    control$iter.max <- 1000
  }
  if("trace" %in% names(control) == FALSE){
    control$trace <- FALSE
  }
  if("constrain" %in% names(control) == FALSE){
    control$constrain <- FALSE
  }

  #### prepare data (scaling, orthogonalization)
  if(control$trace){cat("Scale and center dataset \n")}
  resData <- prepareData.lvm(x, data = data, penalty = x$penalty,
                             regularizationPath = regularizationPath)
  data <- resData$data
  penalty <- resData$penalty
  if(regularizationPath == 1){
    control$start <- resData$start
    x <- resData$x
    names.coef <- resData$names.coef
    orthogonalizer <- resData$orthogonalizer
  }
  
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
  
  if(regularizationPath > 0){
    save_lambda2 <-  penalty$lambda2
  }
  
  control$penalty <- penalty
  
  ## pass parameters for the penalty through the control argument
  control$proxGrad <- list(method = method.proxGrad,
                           step = step,
                           BT.n = BT.n,
                           BT.eta = BT.eta,
                           expX = if(control$constrain){names.coef[x$index$parBelongsTo$cov]}else{NULL},
                           force.descent = force.descent,
                           trace = if(regularizationPath == 0){control$trace}else{FALSE},
                           fixSigma = if(is.null(fixSigma)){FALSE}else{fixSigma}
  )
 
  control$regPath <- list(type = regularizationPath,
                          stepLambda1 = stepLambda1,
                          increasing = increasing,
                          stopParam = stopParam,
                          stopLambda = stopLambda,
                          fixSigma = fixSigma,
                          trace = if(regularizationPath > 0){control$trace}else{FALSE})
  
  ## proximal operator 
  resOperator <- init.proxOperator(lambda1 = penalty$lambda1, 
                                   lambda2 = penalty$lambda2, 
                                   group.penaltyCoef = penalty$group.penaltyCoef,
                                   regularizationPath = regularizationPath)
  control$proxOperator <- resOperator$proxOperator
  control$objectivePenalty <- resOperator$objectivePenalty
  
  #### initialization
  if(is.null(control$start)){
    if(control$trace){cat("Initilization: ")}
    res.init  <- initializer.lvm(x, data = data, names.coef = names.coef, n.coef = n.coef,
                                 penalty = control$penalty, regPath = control$regPath, ...)
  
    if(regularizationPath == 1 || (regularizationPath == 2 && increasing)){
      if(control$trace){cat(" unpenalized LVM \n")}
      control$start <- res.init$lambda0
    }else{
      if(control$trace){cat(" LVM where all penalized coefficient are shrinked to 0 \n")}
      control$start <- res.init$lambdaMax
    }
    
    if(regularizationPath == 2){
      control$regPath$beta_lambdaMax <- res.init$lambdaMax
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
    res$model <- x # restaure the penalized model
  }
  class(res) <- append("plvmfit", class(res))
  
  res$penalty <-  control$penalty[c("names.penaltyCoef", "group.penaltyCoef", "lambda1", "lambda2")]
  res$penalty$regularizationPath <- regularizationPath
  
  #### regularization path
  if(regularizationPath>0){
    res$regularizationPath <- res$opt$message
    res$opt$message <- ""
  }
  
  #### update sigma value
  if( regularizationPath > 0 || fixSigma == TRUE){
    if(control$trace){cat("Estimation of the nuisance parameter ")}
    res <- estimateNuisance.lvm(x, plvmfit = res, data = data, control = control,
                                regularizationPath = regularizationPath)
    if(control$trace){cat("- done \n")}
  }
  
  #### rescale parameter to remove the effect of orthogonalization
  if(regularizationPath == 1){
    setPath(res) <- rescaleCoef_glmPath(Mres = as.matrix(getPath(res)), 
                                        penalty = control$penalty, 
                                        orthogonalizer = orthogonalizer)
  }
  
  #### estimate the best model according to the fit parameter
  if(regularizationPath > 0 && !is.null(fit)){
    cat("Best penalized model according to the",fit,"criteria")
    res <- calcLambda(res, fit = fit, data.fit = data, trace = control$trace)
  }
  
  #### export
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
  
  for(iter_link in penalty$names.penaltyCoef){
     x0 <- rmLink(x0, iter_link)
  }

  ## estimate the model
  suppressWarnings(
    x0.fit <- lava:::estimate.lvm(x = x0, data = data, ...)
  )
  start_lambdaMax[names(coef(x0.fit, level = 9)[,"Estimate"])] <- coef(x0.fit, level = 9)[,"Estimate"]
  
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
#' @param regularizationPath the algorithm used to compute the regularization path. 
#' 
prepareData.lvm <- function(x, data, penalty, regularizationPath){
  
  #### rescale data
  varToScale <- names(which(unlist(lapply(data,is.numeric))))
  if(class(data)[1] != "data.frame"){data <- as.data.frame(data)}
  data[, varToScale] <- scale(data[, varToScale, drop = FALSE])
  # data[, c(varToScale) := lapply(.SD,scale), .SDcols = varToScale]
  
  #### handle factors for the penalty BUT DOES NOT WORK IF SEVERAL APPARTION OF THE SAME VARIABLE
  test.factor <- unlist(lapply(manifest(x), function(name){is.factor(data[[name]])}))
  if(any(test.factor)){ # if some manifest variables are factors
    for(iterFactor in manifest(x)[test.factor]){
      levelTempo <- levels(data[[iterFactor]])
      indexTempo <- grep(paste0("+~",iterFactor,"$"), penalty$names.penaltyCoef)
      varTempo <- penalty$names.penaltyCoef[indexTempo]
      
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
        
        penalty$names.penaltyCoef <- c(setdiff(penalty$names.penaltyCoef,
                                               varTempo),
                                       newTempo) 
        
        penalty$group.penaltyCoef <- c(penalty$group.penaltyCoef[-indexTempo],
                                       rep(penalty$group.penaltyCoef[indexTempo], length(newTempo))
        )
        
      }
    }
  }
  
  #### orthogonalize data
  if(regularizationPath == 1){
    
    resOrtho <- prepareData_glmPath(model = x, data = data, penalty = penalty)
    
    return(list(data = resOrtho$data, 
                penalty = resOrtho$penalty,
                x = resOrtho$model,
                start = resOrtho$mu.X,
                names.coef = coef(resOrtho$model),
                orthogonalizer = resOrtho$orthogonalizer))
    
  }else{
    
    return(list(data = data, 
                penalty = penalty))
    
  }
  
}


#' @title Estimate nuisance parameter
estimateNuisance.lvm <- function(x, plvmfit, data, control, regularizationPath){
  
  control$trace <- FALSE
  names.coefVar <- control$penalty$names.varCoef
    
  if(regularizationPath){
    names.coef <- setdiff(names(getPath(plvmfit, getLambda = NULL, rm.duplicated = TRUE)), names.coefVar)
    fit.coef <- getPath(plvmfit, getLambda = NULL, rm.duplicated = TRUE)[, names.coef, drop = FALSE]
    
  }else{
    names.coef <- setdiff(names( coef(plvmfit)), names.coefVar)
    fit.coef <- setNames(data.frame(rbind(coef(plvmfit)[names.coef, drop = FALSE])), names.coef)
  }
  n.coef <- length(names.coef)
  n.knot <- nrow(fit.coef)
  index.knot <- as.numeric(rownames(fit.coef))
    
  for(iterPath in 1:n.knot){ # build the reduced model associated to each LVM (i.e. fix all penalized parameters to their value)
    xConstrain <- x
    class(xConstrain) <- "lvm"
    
    for(iter_p in 1:n.coef){
      if(names.coef[iter_p] %in% control$penalty$names.penaltyCoef){
        if(fit.coef[iterPath,iter_p]==0){
          xConstrain <- rmLink(xConstrain, var1 = names.coef[iter_p])
        }else{
          xConstrain <- setLink(xConstrain, var1 = names.coef[iter_p], value = fit.coef[iterPath,iter_p])
        }
      }
    }
    
    plvmfit2 <- estimate(xConstrain, data = data, control = control)
    
    if(regularizationPath){
      
      setPath(plvmfit, row = index.knot[iterPath]) <- coef(plvmfit2)[names.coefVar]
      
      if(!is.na(getPath(plvmfit, row = index.knot[iterPath], getLambda = "lambda1.abs", getCoef = NULL))){
        setPath(plvmfit, names = "lambda1", row = index.knot[iterPath]) <- getPath(plvmfit, row = index.knot[iterPath], getLambda = "lambda1.abs", getCoef = NULL) / coef(plvmfit2)[names.coefVar[1]]
      }else if(!is.na(getPath(plvmfit, row = index.knot[iterPath], getLambda = "lambda1", getCoef = NULL))){
        setPath(plvmfit, names = "lambda1.abs", row = index.knot[iterPath]) <- getPath(plvmfit, row = index.knot[iterPath], getLambda = "lambda1", getCoef = NULL) * coef(plvmfit2)[names.coefVar[1]]
      }
      setPath(plvmfit, names = "lambda2", row = index.knot[iterPath]) <- getPath(plvmfit, row = index.knot[iterPath], getLambda = "lambda2.abs", getCoef = NULL) / coef(plvmfit2)[names.coefVar[1]]
      
    }else{
      index1 <- match(names(coef(plvmfit)), names.coefVar, nomatch = 0)>0
      index2 <- match(names(coef(plvmfit2)), names.coefVar, nomatch = 0)>0
      if (!is.null(plvmfit$coef)){ 
        plvmfit$coef[index1,] <- plvmfit2$coef[index2,]
      }
      if (!is.null(plvmfit$opt$estimate)){ 
        plvmfit$opt$estimate[index1] <- plvmfit2$coef[index2,"Estimate"]
      }
      plvmfit$penalty$lambda1.abs <- plvmfit$penalty$lambda1
      plvmfit$penalty$lambda1 <- plvmfit$penalty$lambda1.abs/plvmfit2$coef[which(index2)[1],"Estimate"]
      plvmfit$penalty$lambda2.abs <- plvmfit$penalty$lambda2
      plvmfit$penalty$lambda2 <- plvmfit$penalty$lambda2.abs/plvmfit2$coef[which(index2)[1],"Estimate"]
    }
  }
  
  return(plvmfit)
  
}
