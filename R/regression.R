
"regression<-" <- function(object,...,value) UseMethod("regression<-")
"regression<-.lvm" <- function(object, to=NULL, ..., value) {
  if (!is.null(to)) {    
    regfix(object, to=to, ...) <- value
    return(object)
  } else  {
    if (class(value)[1]=="formula") {

      curvar <- index(object)$var
      
      ##      yx <- all.vars(value)
      lhs <- getoutcome(value)
      X <- attributes(terms(value))$term.labels
      yy <- decomp.specials(lhs)
      res <- lapply(X,decomp.specials)
      xx <- unlist(lapply(res, function(x) x[1]))      

      
      yyf <- lapply(yy,function(y) decomp.specials(y,NULL,"[",fixed=TRUE))
      ys <- unlist(lapply(yyf,function(x) x[1]))      
      object <- addvar(object,ys,...)

      
      exo <- c()
      notexo <- c()
      xxf <- lapply(res,function(x) decomp.specials(x,NULL,"[",fixed=TRUE))
      xs <- unlist(lapply(xxf,function(x) x[1]))
      object <- addvar(object,xs,...)

      
      for (i in 1:length(xs)) {        
        xf <- unlist(strsplit(xx[[i]],"[\\[\\]]",perl=TRUE))
        if (length(xf)>1) {
          xpar <- decomp.specials(xf[2],NULL,":")
          val <- ifelse(xpar[1]=="NA",NA,xpar[1])
          valn <- suppressWarnings(as.numeric(val))
          intercept(object,xs[i]) <- ifelse(is.na(valn),val,valn)
          if (length(xpar)>1) {
            val <- ifelse(xpar[2]=="NA",NA,xpar[2])
            valn <- suppressWarnings(as.numeric(val))
            covariance(object,xs[i]) <- ifelse(is.na(valn),val,valn)
          }
          notexo <- c(notexo,xs[i])
        } else { exo <- c(exo,xs[i]) }
      }

      oldexo <- exogenous(object)
      newexo <- setdiff(exo,c(notexo,curvar,ys))
      exogenous(object) <- union(newexo,setdiff(oldexo,notexo))
            
      for (i in 1:length(ys)) {
        y <- ys[i]
        yf <- unlist(strsplit(yy[i],"[\\[\\]]",perl=TRUE))
        if (length(yf)>1) {
          ypar <- decomp.specials(yf[2],NULL,":")
          val <- ifelse(ypar[1]=="NA",NA,ypar[1])
          valn <- suppressWarnings(as.numeric(val))
          intercept(object,y) <- ifelse(is.na(valn),val,valn)
          if (length(ypar)>1) {
            val <- ifelse(ypar[2]=="NA",NA,ypar[2])
            valn <- suppressWarnings(as.numeric(val))
            covariance(object,y) <- ifelse(is.na(valn),val,valn)
          }
        }
        for (j in 1:length(xs)) {        
          if (length(res[[j]])>1) {
            ##            regfix(object, to=y, from=res[[i]][1],...) <- res[[i]][2]
            regfix(object, to=y[1], from=xs[j],...) <- res[[j]][2]
          } else {
            ##            object <- regression(object,to=y,from=res[[i]][1],...)
            object <- regression(object,to=y[1],from=xs[j],...)
          }
        }
      }
      
##      index(object) <- reindex(object)
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
  function(object=lvm(),to,from,fn=NA,silent=lava.options()$silent,
##           newindex=TRUE,
           ...) {
    if (missing(to)) {
      return(regfix(object))
      ####...
    }        
    if (class(to)[1]=="formula") {
      regression(object,silent=silent,...) <- to
##      functional(object,to) <- fn
      object$parpos <- NULL
      return(object)
    }
##    browser()

    sx <- strsplit(from,"@")
    xx <- sapply(sx, FUN=function(i) i[1])
    ps <- sapply(sx, FUN=function(i) i[2])
    sx <- strsplit(xx,"$",fixed=TRUE)
    xs <- sapply(sx, FUN=function(i) i[1])    
    fix <- as.numeric(sapply(sx, FUN=function(i) i[2]))
    allv <- index(object)$vars
    
    object <- addvar(object, c(to,xs), silent=silent)
    for (i in to)
      for (j in xs) {
        ##object <- addvar(object, c(i,xs), debug=debug, silent=silent)
        ##cancel(object) <- c(i,j)
        ##covfix(object,i,j,exo=TRUE) <- "NA"
        ##        Graph(object) <- addEdge(xs,i,Graph(object))
        object$graph <- addEdge(j,i,object$graph)
        ##        Graph(object) <- addEdge(j,i,Graph(object))
        ##        functional(object,xs,i) <- fn
        if (!is.na(fn))
          functional(object,j,i) <- fn
      }
    
    newexo <- setdiff(xs,c(to,allv))
    exo <- exogenous(object)
    if (length(newexo)>0)
      exo <- unique(c(exo,newexo))
    exogenous(object) <- setdiff(exo,to)
    
    if (lava.options()$debug) {
##      Debug(list("x=",x, " y=",y, " est=",fixed),debug)
      print(object$fix)
    } 
    object$fix[xs,to] <- fix
    object$par[xs,to] <- ps
    object$parpos <- NULL
    
##    if (newindex)
    index(object) <- reindex(object)   
    return(object)
  }
