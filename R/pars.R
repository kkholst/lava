##' @export
`pars` <-
  function(x,...) UseMethod("pars")

##' @S3method pars default
pars.default <- function(x,...) {
  if (!is.null(x$opt$est))
    return(x$opt$est)
  if (!is.null(x$coef))
    return(x$coef)
  return(coef(x))
}

###{{{ pars.multigroupfit
## pars.multigroupfit <- function(x,...) {
##   res <- pars.default(x)
##   lapply(ee$model$lvm,coef))
##   coef()
##}
###}}}

###{{{ pars.lvm

##' @S3method pars lvm
pars.lvm <- function(x, A, P, v, e, ...) {
  parres <- A[index(x)$M1==1]
  diagcorfree <- diag(P)[diag(index(x)$P1)==1]
  parres <- c(parres, diagcorfree)

  if (ncol(A)>1)
  for (i in 1:(ncol(index(x)$P1)-1))
    for (j in (i+1):nrow(index(x)$P1)) {
      if (index(x)$P1[j,i]!=0) {
        parres <- c(parres, P[j,i])
      }
    }
  if (length(parres)>0)
  names(parres) <- paste("p",seq_len(length(parres)),sep="")
  if (!missing(v)) {
    parres <- c( v[which(index(x)$v1==1)], parres)
  }
  if (!missing(e)) {
    parres <- c( parres, e[which(index(x)$e1==1)] )
  }
  return(parres)        
}

###}}} pars.lvm

