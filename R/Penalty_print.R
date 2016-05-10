`print.plvm` <- function(x, ...) {
  
  ## normal display
  lava:::print.lvm(x)
  
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
  cat("Penalty: ", penaltyType,"\n",
      "On     : ", paste(x$penalty$names.penaltyCoef, collapse = " "),"\n")
  }else{
    
    test.lasso <- (x$penalty$group.penaltyCoef<1)*(x$penalty$group.penaltyCoef>0)
    if(any(test.lasso==1)){
    cat("Penalty: ", penaltyType,"\n",
        "On     : ", paste(x$penalty$names.penaltyCoef[test.lasso==1], collapse = " "),"\n")
    }
    
    ls.penalty <- tapply(x$penalty$names.penaltyCoef[test.lasso!=1], x$penalty$group.penaltyCoef[test.lasso!=1],list)
    cat("Penalty: Grouped lasso \n")
    lapply(ls.penalty, function(x){cat("On     :",paste(x, collapse = " "),"\n")})
  }
  cat("\n")
  
  ## export
  invisible(x)
}


##' @export
`print.plvmfit` <-
  function(x,level=2,labels=FALSE,...) {
    
    if(x$penalty$regularizationPath == 0){
    
    Mtempo <- CoefMat(x,labels=labels,level=level,...) 
    ncol.M <- ncol(Mtempo)
    if(x$penalty$lambda1>0 || x$penalty$lambda2>0){
      Mtempo <- rbind(Mtempo, "Penalization:" = rep("", ncol.M))
    }
    if(x$penalty$lambda1>0){
      Mtempo <- rbind(Mtempo, "   L1 lambda" = c(x$penalty$lambda1, rep("",ncol.M-1)))
    }
    if(x$penalty$lambda2>0){
      Mtempo <- rbind(Mtempo, "   L2 lambda" = c(x$penalty$lambda2, rep("",ncol.M-1)))
    }
    
    
    print(Mtempo,quote=FALSE,right=TRUE)
    minSV <- attr(vcov(x),"minSV")
    if (!is.null(minSV) && minSV<1e-12) {
      warning("Small singular value: ", format(minSV))
    }
    pseudo <- attr(vcov(x),"pseudo")
    if (!is.null(pseudo) && pseudo) warning("Singular covariance matrix. Pseudo-inverse used.")
    
    
    }else{
      cat("Regularization path: \n")
      print(x$opt$message)
      cat("estimated using ")
      switch(x$penalty$regularizationPath,
             "1" = cat("LARS algorithm \n"),
             "2" = cat("EPSODE algorithm \n"))
    }
    
    
    invisible(x)
  }