coef2.penalized <- function(x, iter_lambda){
  
  if(is.list(x)){
    
    if(!missing(iter_lambda)){
      x <- x[[iter_lambda]]
    }else{
      res <- lapply(x, function(model){
        c(model@lambda1,
          model@lambda2,
          model@unpenalized, 
          model@penalized, 
          model@nuisance$sigma2)
      })
      
      Mres <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
      colnames(Mres) <- c("lambda1","lambda2",names(x[[1]]@unpenalized),names(x[[1]]@penalized),"sigma2")
      return(Mres)
    }
    
  } 
  
  coef <- c(x@unpenalized, 
            x@penalized, 
            x@nuisance$sigma2)
  return(coef)
}


validLVM <- function(x){
  library(penalized)
  
  lambda1 <- x$penalty$lambda1
  if(is.null(x$control$proxGrad$sigmaMax)){lambda1 <- lambda1*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]}
  lambda2 <- x$penalty$lambda2
  if(is.null(x$control$proxGrad$sigmaMax)){lambda2 <- lambda2*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]}
  
  resPenalized <- penalized(as.formula(paste0(endogenous(x),"~.")), 
                            lambda1 = lambda1, 
                            lambda2 = lambda2, 
                            data = x$data$model.frame, trace = FALSE)
  diffCoef <- coef(x) - coef2.penalized(resPenalized)
  
  return(diffCoef)
}