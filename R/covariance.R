"covariance<-" <- function(object,...,value) UseMethod("covariance<-")

"covariance<-.lvm" <- function(object, var1=NULL, var2=NULL, ..., value) {

  if (!is.null(var1)) {
    if (class(var1)[1]=="formula") {
      lhs <- getoutcome(var1)
      xf <- attributes(terms(var1))$term.labels
      xx <- unlist(lapply(xf, function(x) x[1]))
      if (is.null(lhs)) {
        covfix(object,var1,var2,...) <- value
        object$parpos <- NULL
        return(object)
      }
      else {
        yy <- decomp.specials(lhs)
##        xx <- setdiff(all.vars(var1),yy)
        ##        xx <- var2
##        covfix(object,var1=yy,var2=xx,...) <- value
##        object$parpos <- NULL
##        return(object)
      }
    } else {
      yy <- var1; xx <- var2
    }
    covfix(object,var1=yy,var2=xx,...) <- value
    object$parpos <- NULL
    return(object)
  }
  if (class(value)[1]=="formula") {
    lhs <- getoutcome(value)
    if (is.null(lhs)) {
      return(covariance(object,all.vars(value)))
    }
    yy <- decomp.specials(lhs)

    tt <- terms(value, specials=c("f","v"))
    xf <- attributes(terms(tt))$term.labels
    res <- lapply(xf,decomp.specials)
    nx <- length(xf)
    if (nx==1) {
      if(is.null(attr(tt,"specials")$f) | length(res[[1]])<2) {
        if(is.null(attr(tt,"specials")$v) & is.null(attr(tt,"specials")$f))
##          if(is.null(attr(tt,"specials")$v) | is.null(attr(tt,"specials")$f))
                    
          {
          for (i in yy)
            for (j in res[[1]])
              object <- covariance(object, c(i,j))
        } else {
          covfix(object,var1=yy,var2=NULL) <- res[[1]]
        }
      } else {
        covfix(object,var1=yy,var2=res[[1]][1]) <- res[[1]][2]
      }
      object$parpos <- NULL
      return(object)
    }

    xx <- unlist(lapply(res, function(z) z[1]))
    for (y in yy)
      for (i in 1:length(xx)) {
        if (length(res[[i]])>1) {
          covfix(object, var1=y, var2=res[[i]][1]) <- res[[i]][2]
        } else if ((i+1)%in%attr(tt,"specials")$f | (i+1)%in%attr(tt,"specials")$v) {
          covfix(object, var1=y, var2=NULL) <- res[[i]]
        } else {
          object <- covariance(object,c(y,xx[i]))
        }
      }

    object$parpos <- NULL
    return(object)
    ##       yx <- all.vars(value)
    ##       xx <- attributes(terms(value))$term.labels
    ##       yy <- setdiff(yx,xx)
    ##       return(regression(model,to=yy,from=xx,...))
  }
  ##  if (is.null(var1)) {
  else
    covariance(object,value,...)
##  } else  {
##    covfix(object,var1,var2,...) <- value
##  }
}

`covariance` <-
function(object,var=NULL,...) UseMethod("covariance")

`covariance.lvm` <-
function(object,var=NULL,var2,...) {
  if (!is.null(var)) {
    if (class(var)[1]=="formula") {
      covariance(object,...) <- var
      return(object)
    }

    if (!missing(var2)) {
      for (i in 1:length(var)) {
        c1 <- var[i]
        for (j in 1:length(var2)) {
          c2 <- var2[j]
          object <- addvar(object, c(c1,c2), silent=TRUE)
          ##        cancel(object) <- c(c1,c2)
          object$cov[c1,c2] <- object$cov[c2,c1] <- 1
          object$parpos <- NULL
          index(object) <- reindex(object)
        }
      }
    }
    else {
      for (i in 1:length(var)) {
        c1 <- var[i]
        for (j in i:length(var)) {
          c2 <- var[j]
          object <- addvar(object, c(c1,c2), silent=TRUE)
          ##        cancel(object) <- c(c1,c2)
          object$cov[c1,c2] <- object$cov[c2,c1] <- 1
          object$parpos <- NULL
          index(object) <- reindex(object)
        }
      }
    }
    return(object)
  }
  ##   for (c1 in var)
  ##     for (c2 in var) {
  ##       object <- addvar(object, c(c1,c2), silent=TRUE)
  ##       cancel(object) <- c(c1,c2)
  ##       object$cov[c1,c2] <- object$cov[c2,c1] <- 1
  ##       index(object) <- reindex(object)
  ##     }
  ##   return(object)
  ## }
  else
    return(covfix(object))
}

