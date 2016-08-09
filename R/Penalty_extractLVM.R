`coef0` <-
  function(x,...) UseMethod("coef0")

coef0.plvmfit <- function(x, tol = 0, operator = "<=", penalized = FALSE, value = TRUE){
  
  names.coef <- names(coef(x)) 
  
  if(!is.null(tol)){
    coefTempo <- names.coef[do.call(operator, args = list(abs(coef(x)), tol))]
  }else{
    coefTempo <- names.coef
  }
  
  if(penalized){
    coefTempo <- intersect(coefTempo, x$penalty$name.coef)
  }
  
  if(value){
    return(coef(x)[coefTempo])
  }else{
    return(coefTempo)
  }
}
