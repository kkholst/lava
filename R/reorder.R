## Todo
reorder.lvm <- function(x,vars,...) {
  if (!setequal(vars,vars(x))) stop("Supply vector of all variables")
  m <- x
  addvar(m) <- vars  
}
