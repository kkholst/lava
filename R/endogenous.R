`endogenous` <-
function(x,...) UseMethod("endogenous")

`endogenous.lvmfit` <-
function(x,...) {
  endogenous(Model(x),...)
}

`endogenous.lvm` <-
function(x,top=FALSE,latent=FALSE,...) {
  observed <- manifest(x)
  if (latent) observed <- vars(x)
  if (top) {
    M <- as(Graph(x), Class="matrix")
    res <- c()
    for (i in observed)
      if (!any(M[i,]==1))
        res <- c(res, i)
    return(res)
  }
  exo <- exogenous(x)
  return(setdiff(observed,exo))
}


endogenous.list <- function(x,...) {
  endolist <- c()
  for (i in 1:length(x)) {
    ##    exolist <- c(exolist, exogenous(x[[i]]))
    endolist <- c(endolist, endogenous(x[[i]]))
  }
  endolist <- unique(endolist)
  return(endolist)
##  exolist <- unique(exolist)
##  return(exolist[!(exolist%in%endolist)])
}

`endogenous.multigroup` <-
function(x,...) {
  endogenous(Model(x))
}
