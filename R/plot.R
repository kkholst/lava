###{{{ plot.lvm

`plot.lvm` <-
  function(x,all=FALSE,diag=FALSE,cor=TRUE,labels=FALSE,intercept=FALSE,addcolor=TRUE,plain=FALSE,cex,fontsize1=10,debug=FALSE,noplot=FALSE,graph=list(rankdir="BT"),
         attrs=list(graph=graph),
           unexpr=FALSE,
           parameters=TRUE,
           ...) {
  if (all) {
    diag <- cor <- labels <- intercept <- addcolor <- TRUE
  }
  index(x) <- reindex(x)
##  browser()
  if (!require("Rgraphviz")) stop("package Rgraphviz not available")
  g <- finalize(x,diag=diag,cor=cor,addcolor=addcolor,intercept=intercept,plain=plain,cex=cex,fontsize1=fontsize1,unexpr=unexpr)
  if  (labels) {
    AP <- matrices(x,paste("p",seq_len(index(x)$npar),sep=""))
    mylab <- AP$P; mylab[AP$A!="0"] <- AP$A[AP$A!="0"]
    mylab[!is.na(x$par)] <- x$par[!is.na(x$par)]
    mylab[!is.na(x$covpar)] <- x$covpar[!is.na(x$covpar)]
    g <- edgelabels(g, lab=mylab)
  }
  if (noplot)
    return(g)
  if (debug) {
    plot(g)
  } else {
    
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    dots <- list(...)
    dots$attrs <- attrs
    dots$x <- g
    if (is.null(dots$layoutType) & all(index(m)$A==0))
      dots$layoutType <- "circo"
    g <- do.call("layoutGraph", dots)
    res <- tryCatch(renderGraph(g),error=function(e) NULL)
    options(.savedOpt)
  }
  ## if (!is.null(legend)) {
  ##   op <- par(xpd=TRUE)
  ##   legend(legend, c("Exogenous","Endogenous","Latent","Time to event"),
  ##          pt.cex=1.5, pch=15, lty=0, col=cols[1:4], cex=0.8)
  ##   par(op)
  ## }
  invisible(g)
}

###}}} plot.lvm

###{{{ plot.lvmfit

`plot.lvmfit` <-
  function(x,diag=TRUE,cor=TRUE,type,noplot=FALSE,...) {
    if (!require("Rgraphviz")) stop("package Rgraphviz not available")
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    g <- Graph(x)
    newgraph <- FALSE
    if (is.null(g)) {
      newgraph <- TRUE
      Graph(x) <- finalize(Model(x), diag=TRUE, cor=TRUE, ...)
    }
    if(noplot) return(Graph(x))

    ##cat("Setting up graph...\n")
    if (newgraph) {
      if (missing(type))
        type <- "est"
      x <- edgelabels(x, type=type, diag=diag, cor=cor, ...)
    } else {
      if (!missing(type)) {
        x <- edgelabels(x, type=type, diag=diag, cor=cor, ...)
      }
    }
    g <- Graph(x)
    var <- rownames(covariance(Model(x))$rel)
     if (!cor) {
       delta <- 1
      for (r in 1:(nrow(covariance(Model(x))$rel)-delta) ) {
        for (s in (r+delta):ncol(covariance(Model(x))$rel) ) {
          if (covariance(Model(x))$rel[r,s]==1) {
            g <- removeEdge(var[r],var[s], g)
            g <- removeEdge(var[s],var[r], g)
          }
        }
      }
    }
    if (!diag) {
      for (r in 1:(nrow(covariance(Model(x))$rel)) ) {
        if (isAdjacent(g,var[r],var[r]))
          g <- removeEdge(var[r],var[r],g)
      }
    } 
    m <- Model(x); Graph(m) <- g
    g <- plot(m, diag=diag, cor=cor, ...)
    options(.savedOpt)
    invisible(g)    
  }

###}}} plot.lvmfit

###{{{ plot.multigroup
plot.multigroup <- function(x,diag=TRUE,labels=TRUE,...) {
  k <- x$ngroup
  for (i in 1:k)
    plot(x$lvm[[i]],diag=diag,labels=labels, ...)
}
plot.multigroupfit <- function(x,...) {
  plot(Model(x),...)
}
###}}}

