
"regression<-" <- function(object,...,value) UseMethod("regression<-")
"regression<-.lvm" <- function(object, to=NULL, ..., value) {
  if (!is.null(to)) {    
    regfix(object, to=to, ...) <- value
    return(object)
  } else  {
    if (class(value)[1]=="formula") {
      ##      yx <- all.vars(value)
      lhs <- getoutcome(value)
      xf <- attributes(terms(value))$term.labels
      yy <- decomp.specials(lhs)
      res <- lapply(xf,decomp.specials)
      xx <- unlist(lapply(res, function(x) x[1]))
      object <- addvar(object,yy,...)
      for (y in yy)
        for (i in 1:length(res)) {
          if (length(res[[i]])>1)
            regfix(object, to=y, from=res[[i]][1],...) <- res[[i]][2]
          else
            object <- regression(object,to=y,from=res[[i]][1],...)
        }
      object$parpos <- NULL
      return(object)
##       yx <- all.vars(value)
##       xx <- attributes(terms(value))$term.labels
##       yy <- setdiff(yx,xx)
##       return(regression(object,to=yy,from=xx,...))
    }
    if (!is.list(value) | length(value)>2) stop("Value should contain names of outcome (to) and predictors (from)")
    if (all(c("to","from")%in%names(value))) {

      xval <- value$x; yval <- value$y
    } else {
      yval <- value[[1]]; xval <- value[[2]]
    }
##    print(xval); print(yval)
    regression(object, to=yval, from=xval,...)
  }
}

`regression` <-
  function(object,to,from,...) UseMethod("regression")

`regression.lvm` <-
  function(object=lvm(),to,from,fn=NA,debug=FALSE,silent=FALSE,...) {
    if (missing(to)) {
      return(regfix(object))
      ####...
    }        
    if (class(to)[1]=="formula") {
      regression(object,debug=debug,silent=silent,...) <- to
##      functional(object,to) <- fn
      object$parpos <- NULL
      return(object)
    }
    sx <- strsplit(from,"@")
    xx <- sapply(sx, FUN=function(i) i[1])
    ps <- sapply(sx, FUN=function(i) i[2])
    sx <- strsplit(xx,"$",fixed=TRUE)
    xs <- sapply(sx, FUN=function(i) i[1])    
    fix <- as.numeric(sapply(sx, FUN=function(i) i[2]))

    object <- addvar(object, c(to,xs), debug=debug, silent=silent)    
    for (i in to)
      for (j in xs) {
##        object <- addvar(object, c(i,xs), debug=debug, silent=silent)
##        cancel(object) <- c(i,j)
        covfix(object,i,j,exo=TRUE) <- "NA"
##        Graph(object) <- addEdge(xs,i,Graph(object))
        Graph(object) <- addEdge(j,i,Graph(object))
        ##        functional(object,xs,i) <- fn
        functional(object,j,i) <- fn
      }

    if (debug) {
##      Debug(list("x=",x, " y=",y, " est=",fixed),debug)
      print(object$fix)
    } 
    object$fix[xs,to] <- fix
    object$par[xs,to] <- ps
    object$parpos <- NULL
    index(object) <- reindex(object)   
    return(object)
  }
