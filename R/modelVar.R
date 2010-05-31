
###{{{ modelVar
`modelVar` <-
  function(x,p,...) UseMethod("modelVar")

modelVar.lvmfit <- function(x, p=pars(x), ...) modelVar(Model(x),p=p,...)
modelVar.lvm <- function(x,p,data,debug=FALSE,...) {
  pp <- modelPar(x,p)
  res <- moments(x, p=p, data=data,...)
  attr(res, "pars") <- pp$p
  attr(res, "meanpar") <- pp$meanpar
  res
}
###}}} modelVar
