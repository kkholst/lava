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


validLVM <- function(x, GS = NULL){
  library(penalized)
  
  if(is.null(GS)){
    lambda1 <- x$penalty$lambda1
    lambda1 <- lambda1*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]
    lambda2 <- x$penalty$lambda2
    lambda2 <- lambda2*coef(x)[ paste(endogenous(x),endogenous(x),sep = ",")]
    
    GS <- penalized(as.formula(paste0(endogenous(x),"~.")), 
                              lambda1 = lambda1, 
                              lambda2 = lambda2, 
                              data = x$data$model.frame, trace = FALSE)
  }
  
  diffCoef <- coef(x) - coef2.penalized(GS)
  
  return(diffCoef)
}