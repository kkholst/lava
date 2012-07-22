##' @export
"latent<-" <- function(x,...,value) UseMethod("latent<-")

##' @S3method latent<- lvm
"latent<-.lvm" <- function(x, clear=FALSE,..., value) {
  if (class(value)[1]=="formula") {
    return(latent(x,all.vars(value),clear=clear,...))
  }
  latent(x, var=value, clear=clear,...)
}

##' @export
`latent` <-
function(x,...) UseMethod("latent")

##' @S3method latent lvm
`latent.lvm` <-
function(x,var,clear=FALSE,zero=TRUE,silent=lava.options()$silent,...) {
  if (missing(var)) {    
    latentidx <- unlist(x$latent)
    if (length(latentidx)>0)
       return(names(latentidx))
    else
      return(NULL)
  }
  if (clear) {
    x$noderender$shape[var] <- "rectangle"
    x$latent[var] <- NULL
    if (zero) {
      intfix(x,var) <- NA
    }
  } else {
    if (!all(var%in%vars(x))) {
      addvar(x,silent=silent) <- setdiff(var,vars(x))
    }
    x$noderender$shape[var] <- "ellipse"
    x$latent[var] <- TRUE
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

##' @S3method latent lvmfit
`latent.lvmfit` <-
  function(x,clear=FALSE,...) {
    latent(Model(x),...)
  }

##' @S3method latent list
latent.list <- function(x,...) {
  latlist <- c()
  for (i in 1:length(x)) {
    latlist <- c(latlist, latent(x[[i]]))
  }
  latlist <- unique(latlist)
  return(latlist)
}

##' @S3method latent multigroup
`latent.multigroup` <-
function(x,...) {
  latent(Model(x))
}

