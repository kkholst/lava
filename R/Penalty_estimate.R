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
#' @param correctionStep argument for the EPSODE function (see Penalty_EPSODE.R)   [temporary]
#' @param fast argument for the ISTA function (see Penalty_optims.R) 
#' @param step argument for the ISTA function (see Penalty_optims.R)
#' @param BT.n argument for the ISTA function (see Penalty_optims.R)
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
#' lvm.model <- lvm(list(y~X1+X2+X3+X4))
#' df.data <- sim(lvm.model,300)
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
#' ## 
#' p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1)
#' p1lvm.fit
#' 
#' p1lvm.fit_bis <- estimate(plvm.model,  data = df.data, lambda1 = rp1lvm.fit$opt$message[2,"lambda1"])
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
                          regularizationPath = FALSE, stepLambda1 = 50, correctionStep = FALSE,
                          fast = FALSE, step = 1, BT.n = 20, BT.eta = 0.5, trace = FALSE,
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
    penalty$lambda1 <- lambda1
  }
  if(!missing(lambda2)){
    penalty$lambda2 <- lambda2
  }
  penalty$names.varCoef <- names.coef[x$index$parBelongsTo$cov]
  
  if(regularizationPath > 0){
    save_lambda2 <-  penalty$lambda2
  }
  
  control$penalty <- penalty
  
  ## identifiability issue!!
  if("constrain" %in% names(control) == FALSE){
    control$constrain <- FALSE
  }
  if(regularizationPath == 1){fixSigma <- TRUE}
  
  if(fixSigma == TRUE){
    constrain.sigma <- setNames(1-control$constrain, names.coef[x$index$parBelongsTo$cov[1]])
    constrain.sigma_max <- var(data[,strsplit(names(constrain.sigma), split = ",", fixed = TRUE)[[1]][1]])
    
    if(any(penalty$names.penaltyCoef %in% names(constrain))){
      stop("estimate.plvm: wrong specification of the penalty \n",
           "cannot penalize the variance parameter when fixSigma = TRUE \n",
           "requested penalisation on variance parameters: ",penalty$names.penaltyCoef[penalty$names.penaltyCoef %in% names(constrain)],"\n")
    }
    
  }else{
    
    constrain.sigma <- NULL
    constrain.sigma_max <- NULL
    
  }
  
  ## pass parameters for the penalty through the control argument
  
  control$proxGrad <- list(fast = fast,
                           step = step,
                           BT.n = BT.n,
                           BT.eta = BT.eta,
                           trace = if(regularizationPath == 0){trace}else{FALSE},
                           fixSigma = constrain.sigma,
                           sigmaMax = constrain.sigma_max
  )
  
  control$regPath <- list(type = regularizationPath,
                          stepLambda1 = stepLambda1,
                          correctionStep = correctionStep,
                          trace = if(regularizationPath > 0){trace}else{FALSE})
  
  ## proximal operator 
  if(any(penalty$group.penaltyCoef>=1)){ # grouped lasso
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){ ## proximal operator
      levels.penalty <- unique(test.penalty1)
      for(iter_group in 1:length(levels.penalty)){
        index_group <- which(test.penalty1==levels.penalty[iter_group])
        x[index_group] <- proxE2(x = x[index_group], step = step, lambda = lambda1[index_group])
      }
      return(x)
    }
  }else if(all(penalty$lambda1 == 0) && all(penalty$lambda2 == 0) && regularizationPath == FALSE){ ## no penalty
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){x}
  }else if(all(penalty$lambda2 == 0)){ ## lasso penalty
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){
      mapply(proxL1, x = x, step = step, lambda = lambda1, test.penalty = test.penalty1)
    }
  }else if(all(penalty$lambda1 == 0) &&  regularizationPath == FALSE){ ## ridge penalty
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){
      mapply(proxL2, x = x, step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
  }else{
    proxOperator <- function(x, step, lambda1, lambda2, test.penalty1, test.penalty2){## elastic net penalty
      mapply(proxL2, 
             x = mapply(proxL1, x = x, step = step,  lambda = lambda1, test.penalty = test.penalty1),
             step = step, lambda = lambda2, test.penalty = test.penalty2)
    }
  }
  control$proxOperator <- proxOperator
  
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
    control$start <- initializer.lvm(x, data = data, names.coef = names.coef, n.coef = n.coef,
                                     penalty = control$penalty, regPath = control$regPath, ...)
  }

  #### main
  res <- lava:::estimate.lvm(x = x, data = data, estimator = "penalized", 
                             method = if(regularizationPath == 0){"optim.proxGrad"}else{"optim.regPath"}, 
                             control = control, ...)
  
  #### rescale parameter to remove the effect of orthogonalization
  if(regularizationPath == 1){
    res$opt$message <- rbind(c(0,0,0,coef(do.call(lava:::estimate.lvm, args = c(list(x = x, data = data))))),
                             res$opt$message
    )
    
    res$opt$message <- rescaleRes(Mres = as.matrix(res$opt$message), 
                                  penalty = control$penalty, 
                                  orthogonalizer = resOrtho$orthogonalizer)
    
    
  }
  
  #### export
  res$penalty <-  control$penalty[c("names.penaltyCoef", "group.penaltyCoef", "lambda1", "lambda2")]
  res$penalty$regularizationPath <- regularizationPath
  if(regularizationPath > 0){
    res$opt$message$lambda2 <- save_lambda2
  }
  
  class(res) <- append("plvmfit", class(res))
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
  
  if(regPath$type == 1 || (regPath$type == 2 && regPath$stepLambda1 > 0)){
    initLVM <- try(lava:::estimate.lvm(x = x, data = data, ...))
  }else{
    initLVM <- 0
    class(initLVM) <- "try-error"
  }
  
  if("try-error" %in% class(initLVM) == FALSE){ ## normal model
      
      start <- coef(initLVM)
      
  }else{ ## hight dimensional model
  
  x0 <- x
  start <- setNames(rep(0,n.coef),names.coef)
  
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
  start[names(coef(x0.fit))] <- coef(x0.fit)
  
  }
  
  return(start)
}


#' @title Prepare the data for the glmPath algorithm
prepareDataPath.lvm <- function(model, data, penalty, label = "_pen"){
  
  outcomes <- model$index$endogenous
  n.outcomes <- length(outcomes)
  exogeneous <- model$exogenous
  names.coef <- coef(model)
  n.coef <- length(names.coef)
  
  #### rename penalized variables 
  for(iter_Y in 1:n.outcomes ){
    Y_tempo <- outcomes[iter_Y]
    
    index_coef <- grep(Y_tempo, x = penalty$names.penaltyCoef, fixed = TRUE) # any penalize coefficient related to outcome Y_tempo
    if( length(index_coef)>0 ){
      ls.tempo <- lapply(exogeneous, grep, penalty$names.penaltyCoef[index_coef], fixed = TRUE) # which exogeneous varaible is related to outcome Y_tempo
      names.varData <- exogeneous[unlist(lapply(ls.tempo, function(t){length(t)>0}))] 
      
      data <- cbind(data,
                    setNames(data[, names.varData, drop = FALSE], paste0(names.varData, label, Y_tempo))
      )
      
      # remove link corresponding to a penalized coefficient
      for(iter_link in unlist(ls.tempo)){
        cancel(model) <- as.formula(penalty$names.penaltyCoef[index_coef][iter_link])
      }
      
      # add new link corresponding to the penalized coeffcient but renaæed
      for(iter_link in unlist(ls.tempo)){
        regression(model) <- as.formula(paste0(penalty$names.penaltyCoef[index_coef][iter_link], label, Y_tempo))
      }
      
      # update penalty
      penalty$names.penaltyCoef[index_coef][unlist(ls.tempo)] <- paste0(penalty$names.penaltyCoef[index_coef][unlist(ls.tempo)], label, Y_tempo)
    }
    
  }
  
  #### remove useless variables
  varCoef <- unlist(strsplit(coef(model), split = "~", fixed = TRUE))
  varCoef <- unlist(strsplit(varCoef, split = ",", fixed = TRUE))
  
  indexRM <- which(model$exogenous %in% unique(varCoef) == FALSE)
  if(length(indexRM)>0){
    data <- data[,-which(names(data) %in%  model$exogenous[indexRM])]
    kill(model) <- model$exogenous[indexRM]
  }
  
  #### orthogonalize the dataset
  names.coef <- coef(model)
  lambda1 <- setNames(rep(0, n.coef),names.coef)
  lambda2 <- setNames(rep(0, n.coef),names.coef)
  mu.X <- setNames(rep(0, n.coef),names.coef)
  sd.X <- setNames(rep(1, n.coef),names.coef)
  orthogonalizer <- setNames(vector("list", n.outcomes), outcomes)
  
  for(iter_Y in 1:n.outcomes ){
    Y_tempo <- outcomes[iter_Y]
    
    index_coef <- grep(Y_tempo, x = penalty$names.penaltyCoef, fixed = TRUE) # any penalize coefficient related to outcome Y_tempo
    if( length(index_coef)>0 ){
      
      # orthogonalize data relative to non penalize coefficients
      resOrtho <- orthoData.lvm(model, name.Y = Y_tempo,
                                allCoef = names.coef[grep(Y_tempo, names.coef, fixed = TRUE)], 
                                penaltyCoef = penalty$names.penaltyCoef[index_coef], 
                                data = data)
      
      # update results
      data[, colnames(resOrtho$orthogonalizer)] <- resOrtho$data[, colnames(resOrtho$orthogonalizer),drop = FALSE]
      lambda1[names(resOrtho$lambda1)] <- resOrtho$lambda1
      lambda2[names(resOrtho$lambda2)] <- resOrtho$lambda2
      mu.X[names(resOrtho$lambda1)] <- resOrtho$mu.X
      sd.X[names(resOrtho$lambda1)] <- resOrtho$sd.X
      orthogonalizer[[iter_Y]] <- resOrtho$orthogonalizer
    }
    
  }
  
  #### export
  penalty$sd.X <- sd.X
  penalty$lambda1 <- lambda1
  penalty$lambda2 <- lambda2
  
  return(list(model = model,
              penalty = penalty,
              data = data,
              orthogonalizer = orthogonalizer,
              mu.X = mu.X))
  
}

#' @title Orthogonalize the non-penalized variables relatively to the penalized variables for a regression model
orthoData.lvm <- function(model, name.Y, allCoef, penaltyCoef, data){
  
  ## function
  extractVar <- function(names){
    names.formula <- grep("~", names ,fixed = TRUE)  
    if(length(names.formula)>0){
      names[names.formula] <- unlist(lapply(strsplit(names[names.formula], split = "~", fixed = TRUE),"[",2))
    }
    
    return(names)
  }
  
  ## preparation
  n <- nrow(data)
  n.coef <- length(allCoef)
  names.interceptCoef <- intersect(allCoef,coef(model)[model$index$parBelongsTo$mean])
  n.interceptCoef <- length(names.interceptCoef)
  names.covCoef <- intersect(allCoef,coef(model)[model$index$parBelongsTo$cov])
  if(length(model$latent)>0){
    names.latentCoef <- allCoef[sapply(names(model$latent), grep, x = allCoef, fixed = TRUE)]
  }else{
    names.latentCoef <- NULL
  }
  
  var.penalized <-  extractVar(penaltyCoef)
  var.unpenalized <- setdiff(extractVar(setdiff(allCoef,c(names.covCoef,names.latentCoef))), 
                             c(var.penalized))
  
  ## rebuild data
  X_tempo <- data[,setdiff(c(var.penalized,var.unpenalized),names.interceptCoef), drop = FALSE]
  
  if(n.interceptCoef>0){
    X_tempo <- cbind(matrix(1, nrow = n, ncol = n.interceptCoef),
                     X_tempo)
  }
  names(X_tempo)[1:n.interceptCoef] <- names.interceptCoef
  
  ## distinguish penalized from non penalized
  penalized <-  as.matrix(X_tempo[,var.penalized, drop = FALSE])
  unpenalized <-  as.matrix(X_tempo[,var.unpenalized, drop = FALSE])
  
  ## orthogonlize
  orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
  penalized <- penalized - unpenalized %*% orthogonalizer
  
  ## starting coefficients
  mu.X <- setNames(rep(0,n.coef), allCoef)
  lm.fitted <- lm.fit(y = data[[name.Y]], x = unpenalized)
  mu.X[names(mu.X) %in% penaltyCoef == FALSE] <- c(coef(lm.fitted), var(lm.fitted$residuals)) ### issue with the latent variable here!!
  
  ## scale
  sd.X <- setNames(rep(1,n.coef), allCoef)
  
  varNI.penalized <- setdiff(var.penalized, names.interceptCoef)
  index.penalized <- setdiff(which(names(mu.X) %in% penaltyCoef),
                             which(names(mu.X) %in% names.interceptCoef) )
  
  if(length(varNI.penalized)>0){
    sd.X[index.penalized] <- sqrt(apply(penalized[,varNI.penalized, drop = FALSE], 2, var)*(n-1)/n)
    penalized <- sweep(penalized[,varNI.penalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.penalized])
  }
  
  varNI.unpenalized <- setdiff(var.unpenalized, names.interceptCoef)
  index.unpenalized <- setdiff(which(names(mu.X) %in% penaltyCoef == FALSE), 
                               which(names(mu.X) %in% c(names.interceptCoef, names.covCoef))
  )
  
  if(length(varNI.unpenalized)>0){
    sd.X[index.unpenalized] <- sqrt(apply(unpenalized[,varNI.unpenalized, drop = FALSE], 2, var)*(n-1)/n)
    unpenalized <- sweep(unpenalized[,varNI.unpenalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.unpenalized])
  }
  mu.X <- mu.X * sd.X
  
  ## update initial dataset
  data[,var.penalized] <- penalized
  if(length(setdiff(var.unpenalized, names.interceptCoef))>0){
    data[,setdiff(var.unpenalized, names.interceptCoef)] <- unpenalized
  }
  
  ## lambda
  lambda1 <- setNames(rep(0,n.coef), allCoef)
  lambda1[which(names(mu.X) %in% penaltyCoef)] <- 1/sd.X[which(names(mu.X) %in% penaltyCoef)]
  lambda2 <- setNames(rep(0,n.coef), allCoef)
  lambda2[which(names(mu.X) %in% penaltyCoef)] <- 1/(sd.X[which(names(mu.X) %in% penaltyCoef)]^2)
  
  ## export
  return(list(data = data,
              orthogonalizer = orthogonalizer,
              sd.X = sd.X,
              mu.X = mu.X,
              lambda1 = lambda1,
              lambda2 = lambda2))
}

#' @title Cancel the effect of the orthogonalization on the estimated parameters

rescaleRes <- function(Mres, penalty, orthogonalizer){
  
 
  #### rescale
  names.rescale <- names(penalty$sd.X)
  Mres[,names.rescale] <- sweep(Mres[,names.rescale], MARGIN = 2, FUN = "/", STATS = penalty$sd.X)
  
  #### de-orthogonalize
  name.Y <- names(orthogonalizer)
  covCoef <- penalty$names.varCoef
  n.Y <- length(name.Y)
  
  for(iter_Y in 1:n.Y){
    Y_tempo <- name.Y[iter_Y]
    names.allCoef <- colnames(Mres)[grep(Y_tempo, x = colnames(Mres), fixed = TRUE)] # coefficients related to outcome Y_tempo
    names.penalizedCoef <- intersect(penalty$names.penaltyCoef, names.allCoef) # penalized coefficient related to outcome Y_tempo
    names.unpenalizedCoef <- setdiff(names.allCoef, c(names.penalizedCoef,covCoef)) # penalized coefficient related to outcome Y_tempo
    
    Mres[,names.unpenalizedCoef] <- Mres[,names.unpenalizedCoef] - Mres[,names.penalizedCoef] %*% t(orthogonalizer[[iter_Y]])
  }
  
  return(as.data.frame(Mres))
}