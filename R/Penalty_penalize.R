#### 1- S3 Methods  - init ####
##' @aliases penalize penalize<- penalize.lvm penalize.plvm
##' @param x
##' @param intercept should the intercept be penalized
##' @param regression should the mean parameters be penalized
##' @param variance should the variance parameters be penalized (not possible now)
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

lvm2plvm <- function(x){
  
  x$penalty <- list(names.penaltyCoef = NULL,
                    group.penaltyCoef = NULL,
                    lambda1 = 0, 
                    lambda2 = 0,
                    adaptive = FALSE,
                    V = NULL)
  
  x$penaltyNuclear <- list(FCTobjective = NULL,
                           FCTgradient = NULL,
                           lambdaN  = 0, 
                           name.coef = NULL,
                           name.Y = NULL,
                           name.X = NULL,
                           nrow = NULL,
                           ncol = NULL)
  
  class(x) <- append("plvm", class(x))
  return(x)
}


#### 2- Penalty functions ####

#### Lasso, group lasso, nuclear norm ####
##' @export
`penalize` <-
  function(x,...) UseMethod("penalize")
##' @export
"penalize<-" <- function (x, ..., value) {
  UseMethod("penalize<-", x)
}


`penalize.lvm` <- function(x, value = NULL, ...){
  
  penalize(x, ...) <- value

  return(x)
}

`penalize.plvm` <- `penalize.lvm`

`penalize<-.lvm` <- function(x, ..., value){
  
  ## convert to plvm
  x <- lvm2plvm(x)
  ## main (call `penalty<-.plvm`)
  penalize(x, ...) <- value
  # `penalize<-.plvm`(x, ..., value = value)
  # do.call(`penalize<-.plvm`, args = list(x = x, ..., value = value))
 
  ## export
  return(x)
}

`penalize<-.plvm` <- function(x, intercept = FALSE, regression = TRUE, variance = FALSE, covariance = FALSE, latent = FALSE,
                              lambda1, lambda2, adaptive, V, ..., value){
  
  #### coefficients
  if(!is.null(value)){
    if("formula" %in% class(value)){
      name.X <- all.vars(delete.response(terms(value)))
      name.Y <- setdiff(all.vars(value), name.X)
      value <- paste(name.Y,name.X,sep="~")
    }
    
    if(any(value %in% coef(x) == FALSE)){
      stop("penalty<-.lvm: coefficients to be penalized do not match those of the model\n",
           "unknown coefficients: ",paste(value[value %in% coef(x) == FALSE], collapse = " "),"\n",
           "available coefficients: ",paste(coef(x)[coef(x) %in% value == FALSE], collapse = " "),"\n")
    }
    
    x$penalty$names.penaltyCoef <- value
    
  } else if(is.null(x$penalty$names.penaltyCoef)){

    index.penaltyCoef <- NULL
    if(intercept == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$mean)  
    }
    if(regression == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, x$index$parBelongsTo$reg) 
    }
    if(variance == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, coefVar(x))
    }

    if(covariance == TRUE){
      index.penaltyCoef <- c(index.penaltyCoef, coefCov(x))
    }
    
    ## no penalization on parameters related to the latent variables
    if(latent == FALSE){
      
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
  
   if(!missing(adaptive)){
     x$penalty$adaptive <- as.numeric(adaptive)
   }
   
  #### export
  return(x)
}


#### Nuclear ####
##' @export
`penalizeNuclear` <-
  function(x,...) UseMethod("penalizeNuclear")

##' @export
"penalizeNuclear<-" <- function (x, ..., value) {
  UseMethod("penalizeNuclear<-", x)
}

`penalizeNuclear.lvm` <- function(x, value = NULL, ...){
  
  penalizeNuclear(x, ...) <- value
  
  return(x)
}

`penalizeNuclear.plvm` <- `penalizeNuclear.lvm`

`penalizeNuclear<-.lvm` <- function(x, ..., value){
  
  ## convert to plvm
  x <- lvm2plvm(x)
  
  ## main
  penalizeNuclear(x, ...) <- value
  
  ## export
  return(x)
}

`penalizeNuclear<-.plvm` <- function(x,  coords, objective = NULL, gradient = NULL, lambda = NULL, ..., value){
  
  ## add objective/gradient
  if(is.null(objective) && is.null(x$penaltyNuclear$FCTobjective)){
    x$penaltyNuclear$FCTobjective <- lvGaussian
  }
   
  if(is.null(gradient) && is.null(x$penaltyNuclear$FCTgradient)){
    x$penaltyNuclear$FCTgradient <- scoreGaussian
  }
  
  if(!is.null(objective)){
    x$penaltyNuclear$FCTobjective <- objective
  }
  if(!is.null(gradient)){
    x$penaltyNuclear$FCTgradient <- gradient
  }
  
  ## lambda
  if(!is.null(lambda)){
    x$penaltyNuclear$lambdaN <- lambda  
  }
  
  ## formula
  name.Y <- all.vars(value)[1]
  name.X <- all.vars(value)[-1]
  
  if(name.Y %in% endogenous(x) == FALSE){
    stop("penaltyNuclear: the dependent variable in formula must already be in the model \n")
  }
  if(any(name.X %in% vars(x) == FALSE)){
    stop("penaltyNuclear: the independent variable in formula must be in the model \n",
         "missing name: ",paste(name.X[name.X %in% vars(x) == FALSE],collapse = " "),"\n")
  }
   coords.factor <- apply(coords, 2, function(x){
    as.numeric(as.factor(x))
  })
  ncol <- max(coords.factor)
  nrow <- max(coords.factor)
  
  if(ncol*nrow != NROW(coords)){
    stop("penaltyNuclear: coords does not corresponds to a full matrix \n")
  }
  test.row <- sweep(matrix(coords.factor[,1], nrow = nrow, ncol = ncol),
                    MARGIN = 1, STATS = 1:nrow, FUN = "-")
  if(any(test.row!=0) ){
    stop("penaltyNuclear: coordinates must be ordered \n")
  }
  test.col <- sweep(matrix(coords.factor[,2], nrow = nrow, ncol = ncol),
                    MARGIN = 2, STATS = 1:ncol, FUN = "-")
  if(any(test.col!=0) ){
    stop("penaltyNuclear: coordinates must be ordered \n")
  }
  x$penaltyNuclear$name.coef <- paste0(name.Y,"~",name.X)
  x$penaltyNuclear$name.Y <- name.Y
  x$penaltyNuclear$name.X <- name.X
  x$penaltyNuclear$ncol <- ncol  
  x$penaltyNuclear$nrow <- nrow  
  
  #### export
  return(x)
}

