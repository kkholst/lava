`print.plvm` <- function(x, ...) {
  
  ## normal display
  `print.lvm`(x)
  
  ## additional display
  if(x$penalty$lambda1>0 && x$penalty$lambda2>0){
    penaltyType <- "elastic net"
  }else if(x$penalty$lambda1>0){
    penaltyType <- "L1"
  }else if(x$penalty$lambda1>0){
    penaltyType <- "L2"
  }else{
    penaltyType <- "none"
  }
  cat("Penalty on: ",paste(x$penalty$names.coef, collapse = " "),"\n",
      "Type     : ", penaltyType, "\n")
  cat("\n")
  
  ## export
  invisible(x)
}