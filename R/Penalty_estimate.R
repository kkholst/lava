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
  
  
  # optimizer
  if(proxGrad.method %in% c("ISTA","FISTA") == FALSE){
    stop("estimate.plvm: wrong specification of control$proxGrad.method \n",
         "only  \"ISTA\" and \"FISTA\" mehtods are available \n",
         "control$proxGrad.method: ",proxGrad.method,"\n")
  }  
  
  ## dots  :  cleaning [not completely checked - may be issues here]
  dots <- list(...)
  
  names.dots <- names(dots)
  if(any(names.dots %in% names(penalty))){
    fixedArgs <- setdiff(names(formals("estimate.plvm")),c("..."))
    dotsArgs <- setdiff(names.dots[names.dots %in% names(penalty)], fixedArgs)
    
    if(length(dotsArgs)>0){
      penalty[dotsArgs] <- dots[dotsArgs]
      dots[dotsArgs] <- NULL
    }
  }
  
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
    resOrtho <- orthoData.lvm(x, penalty, df.data)
    data <- resOrtho$data
    dots$control$start <- resOrtho$mu.X
    penalty$sd.X <- resOrtho$sd.X
    penalty$lambda1 <- resOrtho$lambda1
  }
  
  #### initialization
  if(is.null(dots$control$start)){
    x0 <- x
    for(iter_link in penalty$names.penaltyCoef){
      cancel(x0) <- as.formula(iter_link)
    }
    
    x0.fit <- do.call(`estimate.lvm`, args = c(list(x = x0, data = data)))
    dots$control$start <- setNames(rep(0,length(coef(x))),coef(x))
    dots$control$start[names(coef(x0.fit))] <- coef(x0.fit)
  }
  
  #### main
  res <- do.call(`estimate.lvm`, args = c(list(x = x, data = data, estimator = "penalized", 
                                               penalty = penalty, regularizationPath = regularizationPath, 
                                               fix.sigma = fix.sigma, 
                                               method = method, proxGrad.method = proxGrad.method), 
                                          dots)
  )
  
  #### restore the original scale 
  if(regularizationPath){
    res$opt$message <- rbind(res$opt$message,
                             c(0,coef(do.call(`estimate.lvm`, args = c(list(x = x, data = data)))))
    )
    
    M.beta <- as.matrix(res$opt$message)
    index.penalized <- penalty$names.penaltyCoef
    M.beta[,-1] <- sweep(M.beta[,-1], MARGIN = 2, FUN = "/", STATS = resOrtho$sd.X)
    
    index.unpenalized <- setdiff(penalty$names.meanCoef,penalty$names.penaltyCoef)
    M.beta[,index.unpenalized] <- M.beta[,index.unpenalized] - M.beta[,index.penalized] %*% t(resOrtho$orthogonalizer)
    res$opt$message[] <- M.beta
  }
  
  
  
  #### export
  return(res)
}


# beta <- object$beta[unpenalized + seq_len(length(object$beta) - 
#                                             unpenalized)]
# gamma <- object$beta[seq_len(unpenalized)] - drop(orthogonalizer %*% beta)

orthoData.lvm <- function(x, penalty, df.data){
  
  n <- nrow(df.data)
  n.coef <- length(coef(x))
  
  ## function
  extractVar <- function(names){
    names.formula <- grep("~", names ,fixed = TRUE)  
    if(length(names.formula)>0){
      names[names.formula] <- unlist(lapply(strsplit(names[names.formula], split = "~", fixed = TRUE),"[",2))
    }
    
    return(names)
  }
  
  ## rebuild data
  X_tempo <- df.data[,x$index$exogenous, drop = FALSE]
  if(length(penalty$names.interceptCoef)>0){
    X_tempo <- cbind(matrix(1, nrow = nrow(df.data), ncol = length(penalty$names.interceptCoef)),
                     df.data)
  }
  names(X_tempo)[1:length(penalty$names.interceptCoef)] <- penalty$names.interceptCoef
  
  ## distinguish penalized from non penalized
  var.penalized <-  extractVar(penalty$names.penaltyCoef)
  var.unpenalized <-  setdiff(extractVar(penalty$names.meanCoef), var.penalized)
  penalized <-  as.matrix(X_tempo[,var.penalized, drop = FALSE])
  unpenalized <-  as.matrix(X_tempo[,var.unpenalized, drop = FALSE])
  
  ## orthogonlize
  orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
  penalized <- penalized - unpenalized %*% orthogonalizer

  ## center
  mu.X <- setNames(rep(0,n.coef),
                   coef(x))
  lm.fitted <- lm.fit(y = df.data[[x$index$endogenous]], x = unpenalized)
  mu.X[names(mu.X) %in% penalty$names.penaltyCoef == FALSE] <- c(coef(lm.fitted), var(lm.fitted$residuals))
  
  ## scale
  sd.X <- setNames(rep(1,n.coef),
                   coef(x))
  index.penalized <- setdiff(var.penalized, penalty$names.interceptCoef)
  index.penalized2 <- setdiff(which(names(mu.X) %in% penalty$names.penaltyCoef),
                              which(names(mu.X) %in% penalty$names.interceptCoef)
                              )
  if(length(index.penalized)>0){
  sd.X[index.penalized2] <- sqrt(apply(penalized[,index.penalized, drop = FALSE], 2, var)*(n-1)/n)
  penalized <- sweep(penalized[,index.penalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.penalized2])
  }
  
  index.unpenalized <- setdiff(var.unpenalized, penalty$names.interceptCoef)
  index.unpenalized2 <- setdiff(which(names(mu.X) %in% penalty$names.penaltyCoef == FALSE), 
                                which(names(mu.X) %in% c(penalty$names.interceptCoef,penalty$names.varCoef))
  )
  if(length(index.unpenalized)>0){
  sd.X[index.unpenalized2] <- sqrt(apply(unpenalized[,index.unpenalized, drop = FALSE], 2, var)*(n-1)/n)
  unpenalized <- sweep(unpenalized[,index.unpenalized, drop = FALSE], MARGIN = 2, FUN = "/", STATS = sd.X[index.unpenalized2])
  }
  mu.X <- mu.X * sd.X
  
  ## update initial dataset
  df.data[,var.penalized] <- penalized
  if(length(setdiff(var.unpenalized, penalty$names.interceptCoef))>0){
  df.data[,setdiff(var.unpenalized, penalty$names.interceptCoef)] <- unpenalized
  }
  
  ## lambda
  lambda1 <- setNames(rep(0,n.coef),
                      coef(x))
  lambda1[which(names(mu.X) %in% penalty$names.penaltyCoef)] <- 1/sd.X[which(names(mu.X) %in% penalty$names.penaltyCoef)]
  lambda2 <- setNames(rep(0,n.coef),
                      coef(x))
  lambda2[which(names(mu.X) %in% penalty$names.penaltyCoef)] <- 1/(sd.X[which(names(mu.X) %in% penalty$names.penaltyCoef)]^2)
  
  ## export
  return(list(data = df.data,
              orthogonalizer = orthogonalizer,
              sd.X = sd.X,
              mu.X = mu.X,
              lambda1 = lambda1,
              lambda2 = lambda2))
}
