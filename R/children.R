"children" <- function(object,...) UseMethod("children")
"parents" <- function(object,...) UseMethod("parents")

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
