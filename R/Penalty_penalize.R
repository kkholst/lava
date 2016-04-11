#### 1- S3 Methods ####
##' @aliases penalize penalize<- penalize.lvm penalize.plvm
##' @param x
##' @param pen.intercept should the intercept be penalized
##' @param pen.exogenous should the mean parameters be penalized
##' @param pen.variance should the variance parameters be penalized (not possible now)
##' @param lambda1 the lasso penalty
##' @param lambda2 the ridge penalty
##' @param fn_penalty user defined penalty function. Arguments coef, lambda1, lambda2.
##' @param gn_penalty first derivative of the user defined penalty function. Arguments coef, lambda1, lambda2.
##' @param hn_penalty, second derivative of the user defined penalty user defined penalty function. Arguments coef, lambda1, lambda2.
##' @param value
##' @return none
##' @export
##' @examples
##' set.seed(10)
##' n <- 500
##' formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
##' lvm.modelSim <- lvm()
##' regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
##' distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
##' df.data <- sim(lvm.modelSim,n)
##' 
##' lvm.model <- lvm(formula.lvm)
##' plvm.model <- penalize(lvm.model)
##' 
##' #### regularization
##' elvm.L2 <- estimate(plvm.model,  data = df.data, lambda2 = 50)
##' elvm.L1 <- estimate(plvm.model,  data = df.data, lambda1 = 65)
##' elvm.L12 <- estimate(plvm.model,  data = df.data, lambda1 = 50, lambda2 = 50)
##' 
##' #### path regularization
##' elvm.FixedPath <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, fix.sigma = TRUE,
##'                            control = list(data = df.data))
##'                            
##' elvm.L1Fixedpath <- estimate(plvm.model,  data = df.data, lambda1 = elvm.FixedPath$opt$message[3,"lambda"],
##'                         fix.sigma = TRUE)                           
##' coef(L1Fixedpath) - elvm.FixedPath$opt$message[3,-1]
##'
##' elvm.FreePath <- estimate(plvm.model,  data = df.data, regularizationPath = TRUE, control = list(data = df.data))                                                        
##' elvm.L1FreePath <- estimate(plvm.model,  data = df.data, lambda1 = elvm.FreePath$opt$message[3,"lambda"])                           
##' coef(elvm.L1FreePath) - elvm.FreePath$opt$message[3,-1]
##' 
##' @export
`penalize` <-
  function(x,...) UseMethod("penalize")

##' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}


#### 2- Penalty functions ####
`penalize.lvm` <- function(x, ...){
  
  dots <- list(...)
  if("value" %in% names(dots)){
    value <- dots$value
    dots$value <- NULL
    penalize(x, dots) <- dots$value
  }else{
    penalize(x, ...) <- NULL
  }
  
  return(x)
}
`penalize.plvm` <- `penalize.lvm`

`penalize<-.lvm` <- function(x, ..., value){
  
  ## add slot to object
  if("plvm" %in% class(x) == FALSE){
    x$penalty <- list(fn_penalty = NULL,
                      gn_penalty = NULL,
                      hn_penalty = NULL,
                      names.meanCoef = NULL,
                      names.varCoef = NULL,
                      names.interceptCoef = NULL,
                      names.penaltyCoef = NULL,
                      lambda1 = 0, 
                      lambda2 = 0)
    class(x) <- append("plvm", class(x))
  }
  
  ## main (call `penalty<-.plvm`)
  penalize(x, ...) <- value
  # `penalize<-.plvm`(x, ..., value = value)
  # do.call(`penalize<-.plvm`, args = list(x = x, ..., value = value))
 
  ## export
  return(x)
}

`penalize<-.plvm` <- function(x, pen.intercept = FALSE, pen.exogenous = TRUE, pen.variance = FALSE,
                              lambda1, lambda2, fn_penalty, gn_penalty, hn_penalty, ..., value){
  
  index.varCoef <- grep(",", coef(x))
  index.meanCoef <- setdiff(1:length(coef(x)), index.varCoef)
  index.interceptCoef <- setdiff(index.meanCoef, grep("~", coef(x)))
  index.latentCoef <- intersect(grep(names(x$latent), coef(x)), index.meanCoef)
  
  x$penalty$names.varCoef <- coef(x)[index.varCoef]
  x$penalty$names.meanCoef <- coef(x)[index.meanCoef]
  x$penalty$names.interceptCoef <- coef(x)[index.interceptCoef]
  
  ## coefficients
  if(!is.null(value)){
    
    if(any(value %in% coef(x) == FALSE)){
      stop("penalty<-.lvm: coefficients to be penalized do not match those of the model\n",
           "unknown coefficients: ",paste(value[value %in% coef(x) == FALSE], collapse = " "),"\n",
           "available coefficients: ",paste(coef(x)[coef(x) %in% value == FALSE], collapse = " "),"\n")
    }
    
    index.penaltyCoef <- which(coef(x) %in% value)
    x$penalty$names.penaltyCoef <- coef(x)[index.penaltyCoef]
    
  } else if(is.null(x$penalty$names.penaltyCoef)){
  
    index.penaltyCoef <- NULL
    if(pen.intercept == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, index.interceptCoef)
    }
    if(pen.exogenous == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, setdiff(index.meanCoef,c(index.interceptCoef,index.latentCoef)))
    }
    if(pen.variance == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, index.varCoef)
    }
    
    x$penalty$names.penaltyCoef <- coef(x)[index.penaltyCoef]
  } 
  
  ## penalization parameters
  if(!missing(lambda1)){
    x$penalty$lambda1 <- lambda1
  }
  
  if(!missing(lambda1)){
    x$penalty$lambda2 <- lambda2
  }
  
  ## objective function
  if(!missing(fn_penalty)){
    
    x$penalty$fn_penalty <- fn_penalty
    
  } else if(is.null(x$penalty$fn_penalty)){
    
    x$penalty$fn_penalty <-  fn_penalty <- function(coef, lambda1, lambda2){
      fn1 <- lambda1 * sum( abs(coef) )
      fn2 <- lambda2/2 * sum( coef^2 )
      
      return( fn1 + fn2 )
    }
    
  }
  
  ## gradient function
  if(!missing(gn_penalty)){
    
    x$penalty$gn_penalty <- gn_penalty
    
  } else if(is.null(x$penalty$gn_penalty)){
    
    x$penalty$gn_penalty <- function(coef, lambda1, lambda2){
      
      gn1 <- lambda1 * sign(coef)
      gn2 <- lambda2 * abs(coef)
      
      return( gn1 + gn2 )
    }
    
  }
  
  ## hessian function
  if(!missing(hn_penalty)){
    
    x$penalty$hn_penalty <- hn_penalty
    
  } else if(is.null(x$penalty$hn_penalty)){
    
    x$penalty$hn_penalty <- hn_penalty <- function(coef, lambda1, lambda2){
      
      hn1 <- 0
      hn2 <- lambda2 * sign(coef)
      
      return( hn1 + hn2 )
    }
    
  }
  
  ## dots
  dots <- list(...)
  names.dots <- names(dots)
  
  if(length(names.dots) > 0){
    
    
    if(any(names.dots %in% names(x$penalty) == FALSE)){
      fixedArgs <- c("names.coef", "index.coef",
                     setdiff(names(formals("penalty<-.lvm")),c("..."))
      )
      
      stop("penalty<-.lvm: some additional arguments are invalid\n",
           "invalid arguments: ",paste(names.dots[names.dots %in% names(x$penalty) == FALSE], collapse = " "),"\n",
           "valid additional arguments: ", paste(setdiff(names(x$penalty), fixedArgs), collapse = " "),"\n")
      
    }
    x$penalty[names(dots)] <- dots
  }
  
  ## export
  return(x)
}

#### 3- optim functions #### 

penalized_method.lvm <- "proxGrad"#lava:::gaussian_method.lvm # nlminb2

penalized_objective.lvm <- function(x, ...){  # proportional to the log likelihood
  
  obj.UP <- lava:::gaussian_objective.lvm(x,...) # log likelihood
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(!is.null(penalty) && (penalty$lambda1>0 || penalty$lambda2>0) ){ # L1 or L2 penalty
    
    obj.P <- penalty$fn(coef = dots$p[penalty$index.coef],
                        lambda1 = penalty$lambda1,
                        lambda2 = penalty$lambda2)
    
  }else{
    obj.P <- 0
  }
  
  return( obj.UP + obj.P )
}

penalized_gradient.lvm <- function(x, ...){
  
  grad.UP <- gaussian_gradient.lvm(x,...)#lava:::gaussian_gradient.lvm(x,...)
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(!is.null(penalty) && (penalty$lambda1>0 || penalty$lambda2>0) ){ # L1 or L2 penalty
    res_tempo <- penalty$gn(coef = dots$p[penalty$index.coef],
                            lambda1 = penalty$lambda1, 
                            lambda2 = penalty$lambda2)
    
    grad.P <- rep(0, length(dots$p))
    grad.P[penalty$index.coef] <- res_tempo
    
  }else{
    grad.P <- rep(0, length(grad.UP))
  }

  return( grad.UP + grad.P )
  
}

penalized_hessian.lvm <- function(x, ...){ # second order partial derivative
  
  hess.UP <- lava:::gaussian_hessian.lvm(x, ...)
  
  dots <- list(...)
  penalty <- dots$penalty
  
  
  if(!is.null(penalty) && (penalty$lambda1>0 || penalty$lambda2>0) ){ # L1 or L2 penalty
    
    res_tempo <- penalty$hn(coef = dots$p[penalty$index.coef],
                            lambda1 = penalty$lambda1, 
                            lambda2 = penalty$lambda2)
    
    hess.P <- rep(0, length(dots$p))
    hess.P[penalty$index.coef] <- res_tempo
    hess.P <- diag(hess.P)
    
  }else{
    hess.P <- 0
  }
  
  return( hess.UP + hess.P )
  
}

penalized_logLik.lvm <- function(object, ...){ # log likelihood
  
  logLik.UP <- lava:::gaussian_logLik.lvm(object,...)
  
  dots <- list(...)
  penalty <- dots$penalty
  
  if(!is.null(penalty) && (penalty$lambda1>0 || penalty$lambda2>0) ){ # L1 or L2 penalty
    logLik.P <- penalty$fn(coef = dots$p[penalty$names.coef],
                           lambda1 = penalty$lambda1, 
                           lambda2 = penalty$lambda2)
    
  }else{
    logLik.P <- 0
  }
  
  # browser()
  return( logLik.UP + logLik.P )
}

