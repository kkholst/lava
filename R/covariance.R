##' Add covariance structure to Latent Variable Model
##' 
##' Define covariances between residual terms in a \code{lvm}-object.
##' 
##' The \code{covariance} function is used to specify correlation structure
##' between residual terms of a latent variable model, using a formula syntax.
##' 
##' For instance, a multivariate model with three response variables,
##' 
##' \deqn{Y_1 = \mu_1 + \epsilon_1}
##' 
##' \deqn{Y_2 = \mu_2 + \epsilon_2}
##' 
##' \deqn{Y_3 = \mu_3 + \epsilon_3}
##' 
##' can be specified as
##' 
##' \code{m <- lvm(~y1+y2+y3)}
##' 
##' Pr. default the two variables are assumed to be independent. To add a
##' covariance parameter \eqn{r = cov(\epsilon_1,\epsilon_2)}, we execute the
##' following code
##' 
##' \code{covariance(m) <- y1 ~ f(y2,r)}
##' 
##' The special function \code{f} and its second argument could be omitted thus
##' assigning an unique parameter the covariance between \code{y1} and
##' \code{y2}.
##' 
##' Similarily the marginal variance of the two response variables can be fixed
##' to be identical (\eqn{var(Y_i)=v}) via
##' 
##' \code{covariance(m) <- c(y1,y2,y3) ~ f(v)}
##' 
##' To specify a completely unstructured covariance structure, we can call
##' 
##' \code{covariance(m) <- ~y1+y2+y3}
##' 
##' All the parameter values of the linear constraints can be given as the right
##' handside expression of the assigment function \code{covariance<-} if the
##' first (and possibly second) argument is defined as well. E.g:
##' 
##' \code{covariance(m,y1~y1+y2) <- list("a1","b1")}
##' 
##' \code{covariance(m,~y2+y3) <- list("a2",2)}
##' 
##' Defines
##' 
##' \deqn{var(\epsilon_1) = a1}
##' 
##' \deqn{var(\epsilon_2) = a2}
##' 
##' \deqn{var(\epsilon_3) = 2}
##' 
##' \deqn{cov(\epsilon_1,\epsilon_2) = b1}
##' 
##' Parameter constraints can be cleared by fixing the relevant parameters to
##' \code{NA} (see also the \code{regression} method).
##' 
##' The function \code{covariance} (called without additional arguments) can be
##' used to inspect the covariance constraints of a \code{lvm}-object.
##' 
#
##' 
##' @aliases covariance covariance<- covariance.lvm covariance<-.lvm covfix<- covfix covfix<-.lvm covfix.lvm
##' @param object \code{lvm}-object
##' @param var1 Vector of variables names between (or formula)
##' @param var2 Vector of variables names (or formula) defining pairwise
##' covariance between \code{var1} and \code{var2})
##' @param \dots Additional arguments to be passed to the low level functions
##' @param value List of parameter values or (if \code{var1} is unspecified) a
##' @usage
##' \method{covariance}{lvm}(object, var1=NULL, var2=NULL, ...) <- value
##' @return A \code{lvm}-object
##' @author Klaus K. Holst
##' @seealso \code{\link{regression<-}}, \code{\link{intercept<-}},
##' \code{\link{constrain<-}} \code{\link{parameter<-}}, \code{\link{latent<-}},
##' \code{\link{cancel<-}}, \code{\link{kill<-}}
##' @keywords models regression
##' @export
##' @examples
##' 
##' m <- lvm()
##' ### Define covariance between residuals terms of y1 and y2
##' covariance(m) <- y1~y2 
##' covariance(m) <- c(y1,y2)~f(v) ## Same marginal variance
##' covariance(m) ## Examine covariance structure
##'
##' 
`covariance` <- function(object,...) UseMethod("covariance")

##' @export
"covariance<-" <- function(object,...,value) UseMethod("covariance<-")

##' @S3method covariance<- lvm
"covariance<-.lvm" <- function(object, var1=NULL, var2=NULL, ..., value) {

  if (!is.null(var1)) {
    if (class(var1)[1]=="formula") {
      lhs <- getoutcome(var1)
      xf <- attributes(terms(var1))$term.labels
      xx <- unlist(lapply(xf, function(x) x[1]))
      if (length(lhs)==0) {
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
    if (length(lhs)==0) {
      return(covariance(object,all.vars(value),...))
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
              object <- covariance(object, c(i,j),...)
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
          object <- covariance(object,c(y,xx[i]),...)
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

##' @S3method covariance lvm
`covariance.lvm` <-
function(object,var=NULL,var2,exo=FALSE,...) {
  if (!is.null(var)) {
    if (class(var)[1]=="formula") {
      covariance(object,...) <- var
      return(object)
    }
    allvars <- var    
    if (!missing(var2)) {
      if (class(var2)[1]=="formula")
        var2 <- all.vars(var2)
      allvars <- c(allvars,var2)
    }

    xorg <- exogenous(object)
    exoset <- setdiff(xorg,allvars) 
    if (!exo & length(exoset)<length(xorg)) {
##      exogenous(object,mom=TRUE) <- exoset
##      if (length(exoset)==0) exoset <- NA
      exogenous(object) <- exoset
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

