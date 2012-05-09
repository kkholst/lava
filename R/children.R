##' Generic method for extracting children or parent elements of object (e.g. a graph)
##'
##' @title Extract children or parent elements of object
##' @export
##' @aliases children parents
##' @param object Object
##' @param ... Additional arguments
##' @author Klaus K. Holst
"children" <- function(object,...) UseMethod("children")
##' @export
"parents" <- function(object,...) UseMethod("parents")

##' @S3method parents lvmfit
parents.lvmfit <- function(object,...) parents(Model(object),...)

##' @S3method children lvmfit
children.lvmfit <- function(object,...) children(Model(object),...)

##' @S3method parents lvm
parents.lvm <- function(object,var,...) {
  A <- index(object)$A
  if (missing(var)) {
    return(rownames(A))
  }
  if (class(var)[1]=="formula")
    var <- all.vars(var)
  res <- lapply(var, function(v) rownames(A)[A[,v]!=0])
  unique(unlist(res))
}

##' @S3method children lvm
children.lvm <- function(object,var,...) {
  A <- index(object)$A
  if (missing(var)) {
    return(rownames(A))
  }
  if (class(var)[1]=="formula")
    var <- all.vars(var)
  res <- lapply(var, function(v) rownames(A)[A[v,]!=0])
  unique(unlist(res))

}
