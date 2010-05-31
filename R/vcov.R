vcov.lvmfit <- function(object,...) {
  res <- object$vcov
  resnames <- coef(Model(object), mean=object$control$meanstructure)
  colnames(res) <- rownames(res) <- resnames
  return(res)
}
vcov.multigroupfit <- function(object,...) {
  return(object$vcov)
}

