`pars` <-
  function(x,...) UseMethod("pars")

pars.default <- function(x,...) x$coef[,1]
###{{{ pars.multigroupfit
pars.multigroupfit <- function(x,...) {
  x$opt$est
}
###}}}

###{{{ pars.lvm

pars.lvm <- function(x, A, P, v,...) {
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
  names(parres) <- paste("p",1:length(parres),sep="")
  if (!missing(v)) {
    parres <- c( v[which(index(x)$v1==1)], parres)
  }
  return(parres)        
}

###}}} pars.lvm

