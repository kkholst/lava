#' @title Estimate a penalized lvm model
#
#' @param x a penalized lvm model
#' @param data a data.frame containing the data
#' @param lambda1 L1 penalization parameter
#' @param lambda2 L2 penalization parameter
#' @param fn_penalty penalty function relative to the objective function
#' @param gn_penalty penalty function relative to the first derivative of the objective function
#' @param hn_penalty penalty function relative to the second derivative of the objective function
#' @param method optimize to be used
#' @param ... additional arguments to be passed to estimate
#' 
#' @details 
#' Remaining issues with ridge regresion + elastic net to check
#' 
#' @return a
#' 
#' @examples #' 
#' set.seed(10)
#' lvm.model <- lvm(list(y~X1+X2+X3+X4))
#' df.data <- sim(lvm.model,300)
#' plvm.model <- penalize(lvm.model)
#' 
#' ## unpenalized
#' lvm.fit <- estimate(lvm.model, df.data)
#' 
#' ## L1 penalisation
#'  
#' p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1,
#'                control = list(constrain = TRUE, iter.max = 1000))
#' p1lvm.fit
#' 
#' ls.p1lvm <- lapply(c(0,0.1,0.2,0.25,0.5), 
#'                  function(x){estimate(plvm.model,  data = df.data, lambda1 = x, control = list(constrain = TRUE, iter.max = 1000))})
#'                
#' lapply(ls.p1lvm,coef)
#' 
#' ## L2 penalisation
#' 
#' p2lvm.fit <- estimate(plvm.model,  data = df.data, lambda2 = 0.1, 
#'                control = list(constrain = TRUE, iter.max = 1000))
#' p2lvm.fit
#' 
#' ls.p2lvm <- lapply(c(0,0.1,0.2,0.25,0.5), 
#'                  function(x){estimate(plvm.model,  data = df.data, lambda2 = x, control = list(constrain = TRUE, iter.max = 1000))})
#'                
#' lapply(ls.p2lvm,coef)
#' 
#' ## elastic net
#' 
#' p12lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 0.1, lambda2 = 0.1,
#'                control = list(constrain = TRUE, iter.max = 1000))
#' p12lvm.fit
#' 
#' ls.p12lvm <- lapply(c(0,0.1,0.2,0.25,0.5), 
#'                  function(x){estimate(plvm.model,  data = df.data, lambda1 = x, lambda2 = x, control = list(constrain = TRUE, iter.max = 1000))})
#'                
#' lapply(ls.p12lvm,coef)
#' 
#' 
`estimate.plvm` <- function(x, data, lambda1, lambda2, fn_penalty, gn_penalty, hn_penalty, 
                            regularizationPath = FALSE, 
                            fix.sigma = FALSE, 
                            method = penalized_method.lvm, proxGrad.method = "ISTA", ...) {
  
  names.coef <- coef(x)
  n.coef <- length(names.coef)
  
  ## generic penalty
  penalty  <- x$penalty
  
  ## specific values
  if(!missing(lambda1)){
    penalty$lambda1 <- lambda1
  }
  if(!missing(lambda2)){
    penalty$lambda2 <- lambda2
  }
  if(!missing(fn_penalty)){
    penalty$fn_penalty <- fn_penalty
  }
  if(!missing(gn_penalty)){
    penalty$gn_penalty <- gn_penalty
  }
  if(!missing(hn_penalty)){
    penalty$hn_penalty <- hn_penalty
  }
  penalty$names.varCoef <- names.coef[x$index$parBelongsTo$cov]
  
  # optimizer
  if(proxGrad.method %in% c("ISTA","FISTA") == FALSE){
    stop("estimate.plvm: wrong specification of control$proxGrad.method \n",
         "only  \"ISTA\" and \"FISTA\" mehtods are available \n",
         "control$proxGrad.method: ",proxGrad.method,"\n")
  }  
  
  ## dots
  dots <- list(...)
  names.dots <- names(dots)
#   if(any(names.dots %in% names(penalty))){
#     fixedArgs <- setdiff(names(formals("estimate.plvm")),c("..."))
#     dotsArgs <- setdiff(names.dots[names.dots %in% names(penalty)], fixedArgs)
#     
#     if(length(dotsArgs)>0){
#       penalty[dotsArgs] <- dots[dotsArgs]
#       dots[dotsArgs] <- NULL
#     }
#   }
  
  if("control" %in% names.dots == FALSE){
    dots <- list(control = list())
  }
  if("fast" %in% names(dots$control) == FALSE){
    dots$control$fast <- FALSE
  }
  if("iter.max" %in% names(dots$control) == FALSE){
    dots$control$iter.max <- 5000
  }
  
  #### orthogonalization
  if(regularizationPath){
    save_lambda2 <- penalty$lambda2
    save_names.coef <- names.coef
    
    resOrtho <- prepareDataPath.lvm(model = x, data = data, penalty = penalty)

    data <- resOrtho$data
    dots$control$start <- resOrtho$mu.X
    x <- resOrtho$model
    penalty <- resOrtho$penalty
    penalty$lambda2 <- penalty$lambda2 * save_lambda2
    names.coef <- coef(x)
  }
  
  #### initialization
  if(is.null(dots$control$start)){
    x0 <- x
    for(iter_link in penalty$names.penaltyCoef){
      cancel(x0) <- as.formula(iter_link)
    }
    x0.fit <- do.call(`estimate.lvm`, args = c(list(x = x0, data = data)))
    dots$control$start <- setNames(rep(0,n.coef),names.coef)
    dots$control$start[names(coef(x0.fit))] <- coef(x0.fit)
  }
  
  
  #### main
#   regularizationPath <- FALSE
#   penalty$lambda1 <- 100
#   penalty$lambda2 <- 0
  res <- do.call(`estimate.lvm`, args = c(list(x = x, data = data, estimator = "penalized", 
                                               penalty = penalty, regularizationPath = regularizationPath, 
                                               fix.sigma = fix.sigma, 
                                               method = method, proxGrad.method = proxGrad.method), 
                                          dots)
  )

  #### restore the original scale 
  if(regularizationPath){
    res$opt$message$lambda2 <- save_lambda2
    
    res$opt$message <- rbind(res$opt$message,
                             c(0,0,coef(do.call(`estimate.lvm`, args = c(list(x = x, data = data)))))
    )
    
    res$opt$message <- rescaleRes(Mres = as.matrix(res$opt$message), 
                                  penalty = penalty, 
                                  orthogonalizer = resOrtho$orthogonalizer)
    
    names(res$opt$message) <- c("lambda1","lambda2",save_names.coef)
  }
  
  
  
  #### export
  return(res)
}


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
  var.unpenalized <- setdiff(extractVar(allCoef), c(var.penalized,names.covCoef,names.latentCoef))
  
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
  index.penalized <- setdiff(var.penalized, names.interceptCoef)
  index.penalized2 <- setdiff(which(names(mu.X) %in% penaltyCoef),
                              names.interceptCoef)
    
  if(length(index.penalized)>0){
  sd.X[index.penalized2] <- sqrt(apply(penalized[,index.penalized, drop = FALSE], 2, var)*(n-1)/n)
  penalized <- sweep(penalized[,index.penalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.penalized2])
  }
  
  index.unpenalized <- setdiff(var.unpenalized, names.interceptCoef)
  index.unpenalized2 <- setdiff(which(names(mu.X) %in% penaltyCoef == FALSE), 
                                c(names.interceptCoef, model$index$parBelongsTo$cov))
  if(length(index.unpenalized)>0){
  sd.X[index.unpenalized2] <- sqrt(apply(unpenalized[,index.unpenalized, drop = FALSE], 2, var)*(n-1)/n)
  unpenalized <- sweep(unpenalized[,index.unpenalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.unpenalized2])
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

rescaleRes <- function(Mres, penalty, orthogonalizer){
  
  #### rescale
  names.rescale <- names(penalty$sd.X)
  Mres[,names.rescale] <- sweep(Mres[,names.rescale], MARGIN = 2, FUN = "/", STATS = penalty$sd.X)
  
  #### de-orthogonalize
  name.Y <- names(orthogonalizer)
  covCoef <- penalty$names.varCoef
  n.Y <- length(name.Y)
  index.penalized <- penalty$names.penaltyCoef
  
  for(iter_Y in 1:n.Y){
    Y_tempo <- name.Y[iter_Y]
    names.allCoef <- colnames(Mres)[grep(Y_tempo, x = colnames(Mres), fixed = TRUE)] # coefficients related to outcome Y_tempo
    names.penalizedCoef <- intersect(penalty$names.penaltyCoef, names.allCoef) # penalized coefficient related to outcome Y_tempo
    names.unpenalizedCoef <- setdiff(names.allCoef, c(names.penalizedCoef,covCoef)) # penalized coefficient related to outcome Y_tempo
    
    Mres[,names.unpenalizedCoef] <- Mres[,names.unpenalizedCoef] - Mres[,names.penalizedCoef] %*% t(orthogonalizer[[iter_Y]])
  }
  
  return(as.data.frame(Mres))
}