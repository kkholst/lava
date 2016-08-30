##' Generic method for extracting children or parent elements of object (e.g. a graph)
##'
##' @title Extract children or parent elements of object
##' @export
##' @aliases children parents ancestors descendants roots sinks
##' @param object Object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
"children" <- function(object,...) UseMethod("children")
##' @export
"parents" <- function(object,...) UseMethod("parents")
##' @export
"roots" <- function(object,...) UseMethod("roots")
##' @export
"sinks" <- function(object,...) UseMethod("sinks")
##' @export
"descendants" <- function(object,...) UseMethod("descendants")
##' @export
"ancestors" <- function(object,...) UseMethod("ancestors")


##' @export
parents.lvmfit <- function(object,...) parents(Model(object),...)

##' @export
children.lvmfit <- function(object,...) children(Model(object),...)

##' @export
descendants.lvmfit <- function(object,...) descendants(Model(object),...)

##' @export
ancestors.lvmfit <- function(object,...) ancestors(Model(object),...)

##' @export
roots.lvmfit <- function(object,...) roots(Model(object),...)

##' @export
sinks.lvmfit <- function(object,...) sinks(Model(object),...)



##' @export
parents.lvm <- function(object,var,...) {
  A <- index(object)$A
  if (missing(var)) {
    return(rownames(A))
  }
  if (inherits(var,"formula"))
    var <- all.vars(var)
  res <- lapply(var, function(v) rownames(A)[A[,v]!=0])
  unique(unlist(res))
}

##' @export
children.lvm <- function(object,var,...) {
  A <- index(object)$A
  if (missing(var)) {
    return(rownames(A))
  }
  if (inherits(var,"formula"))
    var <- all.vars(var)
  res <- lapply(var, function(v) rownames(A)[A[v,]!=0])
  unique(unlist(res))

}


##' @export
ancestors.lvm <- function(object,x,...) {
   if (inherits(x,"formula")) x <- all.vars(x)
   res <- c()
   left <- setdiff(vars(object),x)
   count <- 0
   child <- x
   while (length(x)>0) {
     count <- count+1
     x <- parents(object,child)
     child <- intersect(x,left)
     res <- union(res,child)
     left <- setdiff(left,child)
   }
   if (length(res)==0) res <- NULL
   return(res)
}

##' @export
descendants.lvm <- function(object,x,...) {
   if (inherits(x,"formula")) x <- all.vars(x)
   res <- c()
   left <- setdiff(vars(object),x)
   count <- 0
   parent <- x
   while (length(x)>0) {
     count <- count+1
     x <- children(object,parent)
     parent <- intersect(x,left)
     res <- union(res,parent)
     left <- setdiff(left,parent)
   }
   if (length(res)==0) res <- NULL
   return(res)
}

##' @export
roots.lvm <- function(object,...) {
    return(exogenous(object,index=FALSE,...))
}

##' @export
sinks.lvm <- function(object,...) {
    return(endogenous(object,top=TRUE,...))
}
