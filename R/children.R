##' Generic method for extracting children or parent elements of object (e.g. a graph)
##'
##' @title Extract children or parent elements of object
##' @export
##' @aliases children parents ancestors descendants
##' @param object Object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
"children" <- function(object,...) UseMethod("children")
##' @export
"parents" <- function(object,...) UseMethod("parents")

##' @export
parents.lvmfit <- function(object,...) parents(Model(object),...)

##' @export
children.lvmfit <- function(object,...) children(Model(object),...)

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
ancestors <- function(m,x,...) {
   if (inherits(x,"formula")) x <- all.vars(x)
   res <- c()
   left <- setdiff(vars(m),x)
   count <- 0
   child <- x
   while (length(x)>0) {
     count <- count+1
     x <- parents(m,child)
     child <- intersect(x,left)
     res <- union(res,child)
     left <- setdiff(left,child)
   }
   if (length(res)==0) res <- NULL
   return(res)
}

##' @export
descendants <- function(m,x,...) {
   if (inherits(x,"formula")) x <- all.vars(x)
   res <- c()
   left <- setdiff(vars(m),x)
   count <- 0
   parent <- x
   while (length(x)>0) {
     count <- count+1
     x <- children(m,parent)
     parent <- intersect(x,left)
     res <- union(res,parent)
     left <- setdiff(left,parent)
   }
   if (length(res)==0) res <- NULL
   return(res)
}

