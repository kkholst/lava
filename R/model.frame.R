model.frame.lvmfit <- function(formula, all=FALSE,...) {
  dots <- list(...)
  mydata <- formula$data$model.frame
  if (all) return(mydata)
  xfix <- colnames(mydata)[(colnames(mydata)%in%parlabels(formula$model0,exo=TRUE))]
  return( mydata[,c(manifest(formula),xfix)] )
}

model.frame.multigroupfit <- function(formula,...) {
  dots <- list(...)
  mydata <- formula$model$data
  return(mydata)
}
