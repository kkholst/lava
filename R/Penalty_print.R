`print.plvm` <- function(x, ...) {
  
  ## normal display
  `print.lvm`(x)
  
  ## additional display
  if(x$penalty$lambda1>0 && x$penalty$lambda2>0){
    penaltyType <- "Elastic net"
  }else if(x$penalty$lambda1>0){
    penaltyType <- "Lasso"
  }else if(x$penalty$lambda1>0){
    penaltyType <- "Ridge"
  }else{
    penaltyType <- "None"
  }
  
  if(all(x$penalty$group.penaltyCoef<1)){
  cat("Penalty on: ",paste(x$penalty$names.penaltyCoef, collapse = " "),"\n",
      "Type      : ", penaltyType, "\n")
  }else{
    test.lasso <- (x$penalty$group.penaltyCoef<1)*(x$penalty$group.penaltyCoef>0)
    if(any(test.lasso==1)){
    x$penalty$group.penaltyCoef[test.lasso==1] <- 0.5
    }
    
    ls.penalty <- tapply(x$penalty$names.penaltyCoef, x$penalty$group.penaltyCoef,list)
    cat("Penalty on: ",paste(ls.penalty[[1]], collapse = " "),"\n")
    lapply(ls.penalty[-1], function(x){cat("          : ",paste(x, collapse = " "),"\n")})
    cat("Type      : Grouped lasso\n")
  }
  cat("\n")
  
  ## export
  invisible(x)
}