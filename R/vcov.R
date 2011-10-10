vcov.lvmfit <- function(object,...) {
  res <- object$vcov
  if ("lvm.missing"%in%class(object)) {
    resnames <- names(coef(object))
  } else {
    resnames <- coef(Model(object), mean=object$control$meanstructure)
  }
  colnames(res) <- rownames(res) <- resnames
  return(res)
}
vcov.multigroupfit <- function(object,...) {
  return(object$vcov)
}

