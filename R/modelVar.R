
###{{{ modelVar

##' @export
`modelVar` <-
  function(x,p,...) UseMethod("modelVar")

##' @S3method modelVar lvmfit
modelVar.lvmfit <- function(x, p=pars(x), ...) modelVar(Model(x),p=p,...)

##' @S3method modelVar lvm
modelVar.lvm <- function(x,p,data,...) {
  pp <- modelPar(x,p)
  res <- moments(x, p=p, data=data,...)
  attr(res, "pars") <- pp$p
  attr(res, "meanpar") <- pp$meanpar
  attr(res, "epar") <- pp$epar
  res
}
###}}} modelVar
