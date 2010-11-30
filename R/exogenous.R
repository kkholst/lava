`exogenous` <-
function(x,...) UseMethod("exogenous")

"exogenous<-" <- function(x,...,value) UseMethod("exogenous<-")

`exogenous<-.lvm` <- function(x,silent=FALSE,
                              xfree=TRUE,
                              ...,value) {
  if (class(value)[1]=="formula") {
    exogenous(x,...) <- all.vars(value)
    return(x)
  }
  not.in <- !(value%in%vars(x))
  if (any(not.in)) {
    addvar(x) <- value[not.in]
  }
  xorg <- exogenous(x)
  x$exogenous <- value
  if (!is.null(value) & xfree) {
    notexo.idx <- which(!(xorg%in%value))
    if (length(notexo.idx)>0) { ##  & mom) {
      if (length(notexo.idx)>1) {
        covariance(x,xorg[notexo.idx],pairwise=TRUE,exo=TRUE) <- NA
      }
      covariance(x,xorg[notexo.idx],vars(x),exo=TRUE) <- NA
      intercept(x,xorg[notexo.idx]) <- NA
    }
  }
##  x$exogenous <- value  
  index(x) <- reindex(x)  
  return(x)
}

`exogenous.lvm` <-
function(x,latent=FALSE,index=TRUE,...) {
  if (!index) {
    if (latent) {
      allvars <- vars(x)
    } else {
      allvars <- manifest(x)
    }
    M <- as(Graph(x), Class="matrix")
    res <- c()
    for (i in allvars)
      if (!any(M[,i]==1) & any(M[i,]==1)) 
        res <- c(res, i)
    return(res)
  }
  if (is.null(x$exogenous)) return(x$exogenous)
  if (all(!is.na(x$exogenous)) & !latent) {
    return(x$exogenous[x$exogenous%in%index(x)$manifest])
  }
  if (!latent)
    return(index(x)$exogenous)
  return(exogenous(x,latent=latent,index=FALSE,...))
}

`exogenous.lvmfit` <-
function(x,...) {
  exogenous(Model(x),...)
}

exogenous.list <- function(x,...) {
  exolist <- c()
  endolist <- c()
  for (i in 1:length(x)) {
    exolist <- c(exolist, exogenous(x[[i]]))
    endolist <- c(endolist, endogenous(x[[i]]))
  }
  endolist <- unique(endolist)
  exolist <- unique(exolist)
  return(exolist[!(exolist%in%endolist)])
}

`exogenous.multigroup` <-
function(x,...) {
  exogenous(Model(x))
}


