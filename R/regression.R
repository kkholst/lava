
##' Add regression association to latent variable model
##'
##' Define regression association between variables in a \code{lvm}-object and
##' define linear constraints between model equations.
##'
##'
##' The \code{regression} function is used to specify linear associations
##' between variables of a latent variable model, and offers formula syntax
##' resembling the model specification of e.g. \code{lm}.
##'
##' For instance, to add the following linear regression model, to the
##' \code{lvm}-object, \code{m}:
##' \deqn{ E(Y|X_1,X_2) = \beta_1 X_1 + \beta_2 X_2}
##' We can write
##'
##' \code{regression(m) <- y ~ x1 + x2}
##'
##' Multivariate models can be specified by successive calls with
##' \code{regression}, but multivariate formulas are also supported, e.g.
##'
##' \code{regression(m) <- c(y1,y2) ~ x1 + x2}
##'
##' defines
##' \deqn{ E(Y_i|X_1,X_2) = \beta_{1i} X_1 + \beta_{2i} X_2 }
##'
##' The special function, \code{f}, can be used in the model specification to
##' specify linear constraints. E.g. to fix \eqn{\beta_1=\beta_2}
##' , we could write
##'
##' \code{regression(m) <- y ~ f(x1,beta) + f(x2,beta)}
##'
##' The second argument of \code{f} can also be a number (e.g. defining an
##' offset) or be set to \code{NA} in order to clear any previously defined
##' linear constraints.
##'
##' Alternatively, a more straight forward notation can be used:
##'
##' \code{regression(m) <- y ~ beta*x1 + beta*x2}
##'
##' All the parameter values of the linear constraints can be given as the right
##' handside expression of the assigment function \code{regression<-} (or
##' \code{regfix<-}) if the first (and possibly second) argument is defined as
##' well. E.g:
##'
##' \code{regression(m,y1~x1+x2) <- list("a1","b1")}
##'
##' defines \eqn{E(Y_1|X_1,X_2) = a1 X_1 + b1 X_2}. The rhs argument can be a
##' mixture of character and numeric values (and NA's to remove constraints).
##'
##' The function \code{regression} (called without additional arguments) can be
##' used to inspect the linear constraints of a \code{lvm}-object.
##'
##' For backward compatibility the "$"-symbol can be used to fix parameters at
##' a given value. E.g. to add a linear relationship between \code{y} and
##' \code{x} with slope 2 to the model \code{m}, we can write
##' \code{regression(m,"y") <- "x$2"}.  Similarily we can use the "@@"-symbol to
##' name parameters. E.g. in a multiple regression we can force the parameters
##' to be equal: \code{regression(m,"y") <- c("x1@@b","x2@@b")}.  Fixed parameters
##' can be reset by fixing (with \$) them to \code{NA}.
##'
##' @aliases regression regression<- regression<-.lvm regression.lvm regfix
##' regfix regfix<- regfix.lvm regfix<-.lvm
##' @param object \code{lvm}-object.
##' @param value A formula specifying the linear constraints or if
##' \code{to=NULL} a \code{list} of parameter values.
##' @param to Character vector of outcome(s) or formula object.
##' @param from Character vector of predictor(s).
##' @param fn Real function defining the functional form of predictors (for
##' simulation only).
##' @param messages Controls which messages are turned on/off (0: all off)
##' @param additive If FALSE and predictor is categorical a non-additive effect is assumed
##' @param y Alias for 'to'
##' @param x Alias for 'from'
##' @param quick Faster implementation without parameter constraints
##' @param \dots Additional arguments to be passed to the low level functions
##' @usage
##' \method{regression}{lvm}(object = lvm(), to, from, fn = NA,
##' messages = lava.options()$messages, additive=TRUE, y, x, value, ...)
##' \method{regression}{lvm}(object, to=NULL, quick=FALSE, ...) <- value
##' @return A \code{lvm}-object
##' @note Variables will be added to the model if not already present.
##' @author Klaus K. Holst
##' @seealso \code{\link{intercept<-}}, \code{\link{covariance<-}},
##' \code{\link{constrain<-}}, \code{\link{parameter<-}},
##' \code{\link{latent<-}}, \code{\link{cancel<-}}, \code{\link{kill<-}}
##' @keywords models regression
##' @examples
##'
##' m <- lvm() ## Initialize empty lvm-object
##' ### E(y1|z,v) = beta1*z + beta2*v
##' regression(m) <- y1 ~ z + v
##' ### E(y2|x,z,v) = beta*x + beta*z + 2*v + beta3*u
##' regression(m) <- y2 ~ f(x,beta) + f(z,beta)  + f(v,2) + u
##' ### Clear restriction on association between y and
##' ### fix slope coefficient of u to beta
##' regression(m, y2 ~ v+u) <- list(NA,"beta")
##'
##' regression(m) ## Examine current linear parameter constraints
##'
##' ## ## A multivariate model, E(yi|x1,x2) = beta[1i]*x1 + beta[2i]*x2:
##' m2 <- lvm(c(y1,y2) ~ x1+x2)
##'
##' @export
"regression<-" <- function(object,...,value) UseMethod("regression<-")

##' @export
regression.formula <- function(object,...) regression(lvm(),object,...)

##' @export
"regression<-.lvm" <- function(object, to=NULL, quick=FALSE, ..., value) {
    dots <- list(...)
    if (length(dots$additive)>0 && !dots$additive && !inherits(value,"formula")) {
        regression(object,beta=value,...) <- to
        return(object)
    }
    if (!is.null(to) || !is.null(dots$y)) {
        regfix(object, to=to, ...) <- value
        return(object)
    } else  {
        if (is.list(value)) {
            for (v in value) {
                regression(object,...) <- v
            }
            return(object)
        }

        if (inherits(value,"formula")) {
            fff <- procformula(object,value,...)
            object <- fff$object
            lhs <- fff$lhs
            xs <- fff$xs
            ys <- fff$ys
            res <- fff$res
            X <- fff$X

            
        if (fff$iscovar) {
            ## return(covariance(object,var1=decomp.specials(lhs[[1]]),var2=X))
            covariance(object) <- toformula(decomp.specials(lhs[[1]]),X)
            return(object)
        }
        if (!is.null(lhs) && nchar(lhs[[1]])>2 && substr(lhs[[1]],1,2)=="v(") {
            v <- update(value,paste(decomp.specials(lhs),"~."))
            covariance(object,...) <- v
            return(object)
        }

        if (length(lhs)==0) {
            index(object) <- reindex(object)
            return(object)
        }

        for (i in seq_len(length(ys))) {
        y <- ys[i]
        for (j in seq_len(length(xs))) {
          if (length(res[[j]])>1) {
            regfix(object, to=y[1], from=xs[j],...) <- res[[j]][2]
          } else {
            object <- regression(object,to=y[1],from=xs[j],...)
          }
        }
      }
      object$parpos <- NULL
      return(object)
    }

    if (!is.list(value) | length(value)>2) stop("Value should contain names of outcome (to) and predictors (from)")
    if (all(c("to","from")%in%names(value))) {

      xval <- value$x; yval <- value$y
    } else {
      yval <- value[[1]]; xval <- value[[2]]
    }
    regression(object, to=yval, from=xval,...)
  }
}

##' @export
`regression` <-
  function(object,to,from,...) UseMethod("regression")

##' @export
`regression.lvm` <-
    function(object=lvm(),to,from,fn=NA,messages=lava.options()$messages,
      additive=TRUE, y,x,value,...) {
        if (!missing(y)) {
            if (inherits(y,"formula")) y <- all.vars(y)
            to <- y
        }
        if (!missing(x)) {
            if (inherits(x,"formula")) x <- all.vars(x)
            from <- x
        }
        if (!additive) {
            if (!inherits(to,"formula")) to <- toformula(to,from)
            x <- attributes(getoutcome(to))$x
            K <- object$attributes$nordinal[x]
            if (is.null(K) || is.na(K)) {
                K <- list(...)$K
                if (is.null(K)) stop("Supply number of categories, K (or use method 'categorical' before calling 'regression').")
                object <- categorical(object,x,...)
            }
            dots <- list(...);
            dots$K <- K
            dots$x <- object
            dots$formula <- to
            dots$regr.only <- TRUE
            object <- do.call("categorical",dots)
            return(object)
        }

        if (missing(to)) {
            return(regfix(object))
        }
        if (inherits(to,"formula")) {
            if (!missing(value)) {
                regression(object,to,messages=messages,...) <- value
            } else {
                regression(object,messages=messages,...) <- to
            }
            object$parpos <- NULL
            return(object)
        }
        if (is.list(to)) {
            for (t in to)
                regression(object,messages=messages,...) <- t
            object$parpos <- NULL
            return(object)
        }

        sx <- strsplit(from,"@")
        xx <- sapply(sx, FUN=function(i) i[1])
        ps <- sapply(sx, FUN=function(i) i[2])
        sx <- strsplit(xx,"$",fixed=TRUE)
        xs <- sapply(sx, FUN=function(i) i[1])
        fix <- char2num(sapply(sx, FUN=function(i) i[2]))
        allv <- index(object)$vars

        object <- addvar(object, c(to,xs), messages=messages, reindex=FALSE)

        for (i in to)
            for (j in xs) {
                object$M[j,i] <- 1
                if (!is.na(fn))
                    functional(object,j,i) <- fn
            }

        if (lava.options()$exogenous) {
            newexo <- setdiff(xs,c(to,allv))
            exo <- exogenous(object)
            if (length(newexo)>0)
                exo <- unique(c(exo,newexo))
            exogenous(object) <- setdiff(exo,to)
        }

        if (lava.options()$debug) {
            print(object$fix)
        }
        object$fix[xs,to] <- fix
        object$par[xs,to] <- ps
        object$parpos <- NULL

        index(object) <- reindex(object)
        return(object)
    }

