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
`penalize.lvm` <- function(x, value = NULL, ...){
  
  penalize(x, ...) <- value

  return(x)
}
`penalize.plvm` <- `penalize.lvm`

`penalize<-.lvm` <- function(x, ..., value){
  
  ## add slot to object
  if("plvm" %in% class(x) == FALSE){
    x$penalty <- list(names.penaltyCoef = NULL,
                      group.penaltyCoef = NULL,
                      lambda1 = 0, 
                      lambda2 = 0,
                      V = NULL)
    class(x) <- append("plvm", class(x))
  }
  
  ## main (call `penalty<-.plvm`)
  penalize(x, ...) <- value
  # `penalize<-.plvm`(x, ..., value = value)
  # do.call(`penalize<-.plvm`, args = list(x = x, ..., value = value))
 
  ## export
  return(x)
}

`penalize<-.plvm` <- function(x, pen.intercept = FALSE, pen.exogenous = TRUE, pen.variance = FALSE, pen.latent = FALSE,
                              lambda1, lambda2, V, ..., value){
  
  #### coefficients
  if(!is.null(value)){
    
    if(any(value %in% coef(x) == FALSE)){
      stop("penalty<-.lvm: coefficients to be penalized do not match those of the model\n",
           "unknown coefficients: ",paste(value[value %in% coef(x) == FALSE], collapse = " "),"\n",
           "available coefficients: ",paste(coef(x)[coef(x) %in% value == FALSE], collapse = " "),"\n")
    }
    
    x$penalty$names.penaltyCoef <- value
    
  } else if(is.null(x$penalty$names.penaltyCoef)){

    index.penaltyCoef <- NULL
    if(pen.intercept == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$mean)  
    }
    if(pen.exogenous == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$reg) 
    }
    if(pen.variance == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$cov)
    }
    
    ## no penalization on parameters related to the latent variables
    if(pen.latent == FALSE){
      
      request <- paste( paste0("^",names(x$latent),"~|~",names(x$latent),"$"), collapse = "|")
      sapply("^u~", grep, x = coef(x), value = TRUE)
      ls.penaltyCoefLatent <- sapply(request, grep, x = coef(x), value = FALSE)
      index.penaltyCoef <- setdiff(index.penaltyCoef, unique(unlist(ls.penaltyCoefLatent)))
    }
      
    x$penalty$names.penaltyCoef <- coef(x)[index.penaltyCoef]
  } 
    #### group penalty if the latent variable is penalized
   names.varLatent <- paste(names(x$latent),names(x$latent),sep = ",")
  if(any(x$penalty$names.penaltyCoef %in% names.varLatent)){
    
    ## check that all paths related to the latent variable are penalized
    for(iter_latent in which(names.varLatent %in% x$penalty$names.penaltyCoef) ){
      index.groupPenalty <- grep(names(x$latent)[iter_latent], coef(x), fixed = TRUE)
      if(any( coef(x)[index.groupPenalty]  %in% x$penalty$names.penaltyCoef == FALSE)){
        message("All paths related to the latent variable ",names(x$latent)[iter_latent]," will be penalized \n",
                "additional penalized parameters: ",paste(setdiff(coef(x)[index.groupPenalty], x$penalty$names.penaltyCoef), collapse = " "),"\n")
        x$penalty$names.penaltyCoef <- union(x$penalty$names.penaltyCoef, coef(x)[index.groupPenalty]) 
      }
    }
    
    ## form groups
    x$penalty$group.penaltyCoef <- seq(0.1, 0.9, length.out = length(x$penalty$names.penaltyCoef))
    for(iter_latent in which(names.varLatent %in% x$penalty$names.penaltyCoef) ){
    index.groupPenalty <- grep(names(x$latent)[iter_latent], x$penalty$names.penaltyCoef, fixed = TRUE)
    x$penalty$group.penaltyCoef[index.groupPenalty] <- iter_latent
    }
    
  }else{
    x$penalty$group.penaltyCoef <- seq(0.1, 0.9, length.out = length(x$penalty$names.penaltyCoef))
  }
   
   #### V matrix
   if(!missing(V)){
     x$penalty$V <- V
   }else if(is.null(x$penalty$V)){
     V <- matrix(0, nrow = length(coef(x)), ncol = length(coef(x)))
     colnames(V) <- coef(x)
     rownames(V) <- coef(x)
     diag(V)[x$penalty$names.penaltyCoef] <- 1
     x$penalty$V <- V
  }
  
   #### penalization parameters
  if(!missing(lambda1)){
    x$penalty$lambda1 <- as.numeric(lambda1)
  }
  
  if(!missing(lambda1)){
    x$penalty$lambda2 <- as.numeric(lambda2)
  }
  
  #### export
  return(x)
}






if(newcoef[1] %in% exogenous(lvmfit)){
  regression(model) <- as.formula(paste(rev(newcoef), collapse = "~"))
}else if(newcoef[2] %in% exogenous(lvmfit)){
  regression(model) <- as.formula(paste(newcoef, collapse = "~"))
}else{
  covariance(model) <- as.formula(paste(newcoef, collapse = "~"))  
}

#### 3- optim functions #### 

penalized_method.lvm <- "proxGrad"#lava:::gaussian_method.lvm # nlminb2

penalized_objective.lvm <- lava:::gaussian_objective.lvm

penalized_objectivePen.lvm <- function(x, lambda1, lambda2){

  obj.P <- sum( lambda1 *  abs(x) + lambda2/2 * x^2 )
  
  return( obj.P )
}


penalized_gradient.lvm <- lava:::gaussian_gradient.lvm

penalized_hessian.lvm <- lava:::gaussian_hessian.lvm
  
penalized_logLik.lvm <- lava:::gaussian_logLik.lvm
