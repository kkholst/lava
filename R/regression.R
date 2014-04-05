##' Add regression association to latent variable model
##' 
##' Define regression association between variables in a \code{lvm}-object and
##' define linear constraints between model equations.
##' 
#
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
##' mixture of character and numeric values (and NA's).
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
##' @param silent Logical variable which indicates whether messages are turned
##' on/off.
##' @param quick Faster implementation without parameter constraints
##' @param \dots Additional arguments to be passed to the low level functions
##' @usage
##' \method{regression}{lvm}(object = lvm(), to, from, fn = NA,
##' silent = lava.options()$silent, ...)
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
##' 	  
##'
##' @export
"regression<-" <- function(object,...,value) UseMethod("regression<-")

##' @S3method regression<- lvm
"regression<-.lvm" <- function(object, to=NULL, quick=FALSE, ..., value) {
  if (!is.null(to)) {    
    regfix(object, to=to, ...) <- value
    return(object)
  } else  {
    if (is.list(value)) {
      for (v in value) {
        regression(object,...) <- v        
      }
      return(object)
    }

    if (class(value)[1]=="formula") {
        yx <- lapply(strsplit(as.character(value),"~"),function(x) gsub(" ","",x))[-1]
        iscovar <- FALSE
        if (length(yx)==1) {        
            lhs <- NULL; xidx <- 1
        } else {
            lhs <- yx[1]; xidx <- 2
            if (yx[[xidx]][1]=="") {
                yx[[xidx]] <- yx[[xidx]][-1]
                iscovar <- TRUE
            }            
        }
        X <- strsplit(yx[[xidx]],"+",fixed=TRUE)[[1]]
        if (iscovar) {
            ## return(covariance(object,var1=decomp.specials(lhs[[1]]),var2=X))
            covariance(object) <- toformula(decomp.specials(lhs[[1]]),X)
            return(object)
        }
        if (!is.null(lhs) && nchar(lhs[[1]])>2 && substr(lhs[[1]],1,2)=="v(") {
            v <- update(value,paste(decomp.specials(lhs),"~."))
            covariance(object,...) <- v
            return(object)
        }        

      curvar <- index(object)$var
      res <- lapply(X,decomp.specials,pattern2="[*]",reverse=TRUE)
      xx <- unlist(lapply(res, function(x) x[1]))      

      notexo <- c()
      if (length(lhs)>0) {
        yy <- decomp.specials(lhs)
        yyf <- lapply(yy,function(y) decomp.specials(y,NULL,pattern2="[",fixed=TRUE))
        ys <- unlist(lapply(yyf,function(x) x[1]))      
        object <- addvar(object,ys,reindex=FALSE,...)
        notexo <- ys
      }
      
      exo <- c()
      xxf <- lapply(as.list(xx),function(x) decomp.specials(x,NULL,pattern2="[",fixed=TRUE))
      xs <- unlist(lapply(xxf,function(x) x[1]))

      ## Remove intercepts?
      rmint <- na.omit(match("-1",xs))
      if (length(rmint)>0) intercept(object,ys) <- 0
      xs <- setdiff(xs,c("-1","1"))

      object <- addvar(object,xs,reindex=FALSE ,...)
      
      for (i in seq_len(length(xs))) {        
        xf <- unlist(strsplit(xx[[i]],"[\\[\\]]",perl=TRUE))
        if (length(xf)>1) {
          xpar <- strsplit(xf[2],":")[[1]]
          if (length(xpar)>1) {
            val <- ifelse(xpar[2]=="NA",NA,xpar[2])
            valn <- suppressWarnings(as.numeric(val))
            covariance(object,xs[i]) <- ifelse(is.na(valn),val,valn)
          }
          val <- ifelse(xpar[1]=="NA",NA,xpar[1])
          valn <- suppressWarnings(as.numeric(val))
          intercept(object,xs[i]) <- ifelse(is.na(valn),val,valn)
          notexo <- c(notexo,xs[i])
        } else { exo <- c(exo,xs[i]) }
      }

      
      if (length(lhs)==0) {
        index(object) <- reindex(object)
        return(object)
      }

      if (lava.options()$exogenous) {
        oldexo <- exogenous(object)
        newexo <- setdiff(exo,c(notexo,curvar,ys))
        exogenous(object) <- union(newexo,setdiff(oldexo,notexo))
      }

      for (i in seq_len(length(ys))) {
        y <- ys[i]
        yf <- unlist(strsplit(yy[i],"[\\[\\]]",perl=TRUE))
        if (length(yf)>1) {
          ypar <- strsplit(yf[2],":")[[1]]
          if (length(ypar)>1) {
            val <- ifelse(ypar[2]=="NA",NA,ypar[2])
            valn <- suppressWarnings(as.numeric(val))
            covariance(object,y) <- ifelse(is.na(valn),val,valn)
          }
          val <- ifelse(ypar[1]=="NA",NA,ypar[1])
          valn <- suppressWarnings(as.numeric(val))
          intercept(object,y) <- ifelse(is.na(valn),val,valn)
        }
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

##' @S3method regression lvm
`regression.lvm` <-
  function(object=lvm(),to,from,fn=NA,silent=lava.options()$silent,
           ...) {
    if (missing(to)) {
      return(regfix(object))
    }        
    if (class(to)[1]=="formula") {
      regression(object,silent=silent,...) <- to
      object$parpos <- NULL
      return(object)
    }
    if (is.list(to)) {
      for (t in to)
        regression(object,silent=silent,...) <- t
      object$parpos <- NULL
      return(object)
    }
    
    sx <- strsplit(from,"@")
    xx <- sapply(sx, FUN=function(i) i[1])
    ps <- sapply(sx, FUN=function(i) i[2])
    sx <- strsplit(xx,"$",fixed=TRUE)
    xs <- sapply(sx, FUN=function(i) i[1])    
    fix <- as.numeric(sapply(sx, FUN=function(i) i[2]))
    allv <- index(object)$vars
    
    object <- addvar(object, c(to,xs), silent=silent,reindex=FALSE)
    
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
