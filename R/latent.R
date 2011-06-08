"latent<-" <- function(x,...,value) UseMethod("latent<-")

"latent<-.lvm" <- function(x, clear=FALSE,..., value) {
  if (class(value)[1]=="formula") {
    return(latent(x,all.vars(value),clear=clear,...))
  }
  latent(x, var=value, clear=clear,...)
}

`latent` <-
function(x,...) UseMethod("latent")

`latent.lvm` <-
function(x,var,clear=FALSE,zero=TRUE,silent=FALSE,...) {
  if (missing(var)) {
    latentidx <- unlist(nodeData(Graph(x), attr="latent"))
    if (length(latentidx)>0)
      return(names(latentidx)[latentidx])
    else
      return(NULL)
  }
  if (clear) {
    x <- addattr(x,attr="shape",var=var,val="rectangle")
    nodeData(Graph(x), var, attr="latent") <- FALSE
    if (zero) {
      intfix(x,var) <- NA
    }
  } else {
    if (!all(var%in%vars(x))) {
      addvar(x,silent=silent) <- setdiff(var,vars(x))
    }
    x <- addattr(x,attr="shape",var=var,val="ellipse")
    nodeData(Graph(x), var, attr="latent") <- TRUE
    if (zero & tolower(lava.options()$param)%in%c("hybrid","absolute")) {
      intercept(x,var) <- 0
    }
  }
  
  xorg <- exogenous(x)
  exoset <- setdiff(xorg,var) 
  if (length(exoset)<length(xorg)) {
    exogenous(x) <- exoset
  }  
  
  index(x) <- reindex(x)
  return(x)
}

`latent.lvmfit` <-
  function(x,clear=FALSE,...) {
    latent(Model(x),...)
  }

latent.list <- function(x,...) {
  latlist <- c()
  for (i in 1:length(x)) {
    latlist <- c(latlist, latent(x[[i]]))
  }
  latlist <- unique(latlist)
  return(latlist)
}

`latent.multigroup` <-
function(x,...) {
  latent(Model(x))
}

