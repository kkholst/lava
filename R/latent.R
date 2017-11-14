##' @export
"latent<-" <- function(x,...,value) UseMethod("latent<-")

##' @export
"latent<-.lvm" <- function(x, clear=FALSE,..., value) {
  if (inherits(value,"formula")) {
    return(latent(x,all.vars(value),clear=clear,...))
  }
  latent(x, var=value, clear=clear,...)
}

##' @export
`latent` <-
function(x,...) UseMethod("latent")

##' @export
`latent.lvm` <- function(x,var,clear=FALSE,messages=lava.options()$messages,...) {
    if (missing(var)) {
        latentidx <- unlist(x$latent)
        if (length(latentidx)>0)
            return(names(latentidx))
        else
            return(NULL)
    }
    if (inherits(var,"formula")) var <- all.vars(var)
    if (clear) {
        x$noderender$shape[var] <- "rectangle"
        x$latent[var] <- NULL
        ## intfix(x,var) <- NA
    } else {
        if (!all(var%in%vars(x))) {
            addvar(x,messages=messages,reindex=FALSE,) <- setdiff(var,vars(x))
        }
        x$noderender$shape[var] <- "ellipse"
        x$latent[var] <- TRUE
        ord <- intersect(var,ordinal(x))
        if (length(ord)>0) ordinal(x,K=NULL) <- ord
    }
    
    xorg <- exogenous(x)
    exoset <- setdiff(xorg,var)
    if (length(exoset)<length(xorg)) {
        exogenous(x) <- exoset
    }
    
    index(x) <- reindex(x)
    return(x)
}

##' @export
`latent.lvmfit` <-
  function(x,clear=FALSE,...) {
    latent(Model(x),...)
  }

##' @export
latent.list <- function(x,...) {
  latlist <- c()
  for (i in seq_along(x)) {
    latlist <- c(latlist, latent(x[[i]]))
  }
  latlist <- unique(latlist)
  return(latlist)
}

##' @export
`latent.multigroup` <-
function(x,...) {
  latent(Model(x))
}
