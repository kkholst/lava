###{{{ plot.lvm

`plot.lvm` <-
  function(x,diag=FALSE,cor=TRUE,labels=FALSE,intercept=FALSE,addcolor=TRUE,plain=FALSE,cex,fontsize1=10,noplot=FALSE,graph=list(rankdir="BT"),
         attrs=list(graph=graph),
           unexpr=FALSE,
           addstyle=TRUE,Rgraphviz=lava.options()$Rgraphviz,
           ...) {
  index(x) <- reindex(x)
  if (length(index(x)$vars)<2) stop("Not available for models with fewer than two variables")
##  browser()
  if (!Rgraphviz || (!require("Rgraphviz"))) {
    if (!require("igraph"))
      stop("package 'Rgraphviz' or 'igraph' not available")
    g <- igraph.lvm(x,...)
    if (noplot) return(g)
    plot(g)
    return(invisible(g))
  } 

  g <- finalize(x,diag=diag,cor=cor,addcolor=addcolor,intercept=intercept,plain=plain,cex=cex,fontsize1=fontsize1,unexpr=unexpr,addstyle=addstyle)
  if  (labels) {
    AP <- matrices(x,paste("p",seq_len(index(x)$npar),sep=""))
    mylab <- AP$P; mylab[AP$A!="0"] <- AP$A[AP$A!="0"]
    mylab[!is.na(x$par)] <- x$par[!is.na(x$par)]
    mylab[!is.na(x$covpar)] <- x$covpar[!is.na(x$covpar)]
    g <- edgelabels(g, lab=mylab)
  }
  if (lava.options()$debug) {
    plot(g)
  } else {

##    graphRenderInfo(g)$recipEdges <- "distinct"
    .savedOpt <- options(warn=-1) ## Temporarily disable warnings as renderGraph comes with a stupid warning when labels are given as "expression"
    dots <- list(...)
    dots$attrs <- attrs
    dots$x <- g
    dots$recipEdges <- "distinct"
    if (is.null(dots$layoutType) & all(index(x)$A==0))
      dots$layoutType <- "circo"
    g <- do.call("layoutGraph", dots)
    if (noplot)
      return(g)    
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

###}}} plot.lvmz

###{{{ plot.lvmfit

`plot.lvmfit` <-
  function(x,diag=TRUE,cor=TRUE,type,noplot=FALSE,...) {
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

###{{{ igraph.lvm

igraph.lvm <- function(x,layout=igraph::layout.kamada.kawai,...) {
  require("igraph")
  oC <- covariance(x)$rel
  for (i in 1:(nrow(oC)-1))
    for (j in (i+1):nrow(oC)) {
      if (oC[i,j]!=0) {
        x <- regression(x,vars(x)[i],vars(x)[j])
        x <- regression(x,vars(x)[j],vars(x)[i])
      }
    }  
  g <- igraph.from.graphNEL(Graph(x))  
  V(g)$color <- "lightblue"
  V(g)$label <- vars(x)
  for (i in match(latent(x),V(g)$name)) {
    V(g)$shape[i] <- "circle"
    V(g)$color[i] <- "green"
  }
  endo <- index(x)$endogenous
  for (i in match(endo,V(g)$name)) {
    V(g)$color[i] <- "orange"
  }
  E(g)$label <- as.list(rep("",length(E(g))))
  oE <- edgelabels(x)
  for (i in 1:length(E(g))) {
    st <- as.character(oE[i])
    if (length(st)>0)
      E(g)$label[[i]] <- st
  }
  g$layout <- layout(g)
  return(g)  
}

###}}} igraph.lvm
