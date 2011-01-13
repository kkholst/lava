
confint.lvmfit <- function(object,parm=1:length(coef(object)),level=0.95,profile=FALSE,...) {
  if (!profile) {
    return(confint.default(object,parm=parm,level=level,...))
  }
  if (is.character(parm)) {
    parm <- parpos(Model(object),p=parm)
    parm <- parm[attributes(parm)$ord]
  } 
  res <- c()
  for (i in parm) {
    res <- rbind(res, profci.lvmfit(object,parm=i,level=level,profile=profile))
  }  
  rownames(res) <- names(coef(object))[parm]
  colnames(res) <- paste((c(0,1) + c(1,-1)*(1-level)/2)*100,"%")
  return(res) 
}


confint.multigroupfit <- function(object,parm=1:length(pars(object)),level=0.95,...) {
  if (is.character(parm)) {
    
  }
  
}
