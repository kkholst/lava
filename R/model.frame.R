model.frame.lvmfit <- function(formula, all=FALSE,...) {
  dots <- list(...)
  mydata <- formula$data$model.frame
  if (all) return(mydata)
  xfix <- colnames(mydata)[(colnames(mydata)%in%parlabels(formula$model0))]
  return( mydata[,c(manifest(formula),xfix)] )
}
