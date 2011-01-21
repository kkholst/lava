
confint.lvmfit <- function(object,parm=1:length(coef(object)),level=0.95,profile=FALSE,...) {
  if (is.character(parm)) {
    parm <- parpos(Model(object),p=parm)
    parm <- parm[attributes(parm)$ord]
  } 
  if (!profile) {
    return(confint.default(object,parm=parm,level=level,...))
  }
  res <- c()
  for (i in parm) {
    res <- rbind(res, profci.lvmfit(object,parm=i,level=level,profile=profile))
  }  
  rownames(res) <- names(coef(object))[parm]
  colnames(res) <- paste((c(0,1) + c(1,-1)*(1-level)/2)*100,"%")
  return(res) 
}


confint.multigroupfit <- function(object,parm=1:length(pars(object)),level=0.95,
                                  estimates=TRUE,...) {
  p <- 1-(1-level)/2
  res <- cbind(pars(object),pars(object)) + qnorm(p)*cbind(-1,1)%x%diag(vcov(object))^0.5
  colnames(res) <- paste(c(1-p,p)*100,"%",sep="")
  rownames(res) <- parpos(object); rownames(res)[is.na(rownames(res))] <- ""
  if (estimates) res <- cbind(coef(object,level=0)[,c(1,2,4)],res)
  res[parm,,drop=FALSE]
}
  
