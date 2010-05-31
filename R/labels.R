###{{{ labels
`labels<-` <- function(object,...,value) UseMethod("labels<-")

`labels<-.default` <- function(object,...,value) {
  labels(object,value)
}

labels.graphNEL <- function(object,lab=NULL,...) {
  if (is.null(lab))
    return(nodeRenderInfo(object)$label)  
  nodeRenderInfo(object) <- list(label=lab)
  names(nodeRenderInfo(object)$label) <- nodes(object);
  return(object)
}

labels.lvmfit <- function(object,lab=NULL,...) {
  if (is.null(lab))
    return(labels(Graph(object)))  
  labels(Graph(object)) <- lab
  return(object)
}

`labels.lvm` <- function(object,lab=NULL,...) {
  gr <- Graph(object)
  if (is.null(lab))
    return(nodeRenderInfo(gr)$label)
##  print(lab)
##  nodeRenderInfo(gr)[] <- list(label=lab)
 ## names(nodeRenderInfo(gr)$label) <- vars(object);
  if (is.null(nodeRenderInfo(gr)$label))
    nodeRenderInfo(gr)$label <- lab
  else 
    nodeRenderInfo(gr)$label[names(lab)] <- lab
  Graph(object) <- gr
  return(object)
}
###}}} labels

###{{{ edgelabels

"edgelabels<-.lvmfit" <- function(object,to,from,est=TRUE,edges=NULL,cex=1,...,value) {
  if (is.null(edges))  {
    if (class(to)[1]=="formula") {
      yy <- getoutcome(to)
      from <- setdiff(all.vars(to),yy)
      to <- yy     
    }
    edges <- paste(from,to,sep="~")
  }
    
  edges. <- paste("\"", edges, "\"", sep="")
  fromto <- edge2pair(edges)
  val <- c()
  for (i in 1:length(edges)) {
    val <- c(val,
             formatC(effects(object,from=fromto[[i]][1],to=fromto[[i]][2],silent=TRUE)$directef[[1]])
             )
  }
  if (est)
    mytext <- paste("c(", paste(paste(edges.,"=expression(",as.character(value),"==\"",val,"\")",sep=""),collapse=","),")")
  else
    mytext <- paste("c(", paste(paste(edges.,"=expression(",as.character(value),")",sep=""),collapse=","),")")
  edgeRenderInfo(Graph(object))$label <- eval(parse(text=mytext))
  edgeRenderInfo(Graph(object))$cex[edges] <- cex  
  return(object)
}

edgelabels.lvmfit <- function(object,value,type,pthres,...) {
  if (!missing(value)) {
    edgelabels(object,...) <- value
    return(object)
  }
  if (missing(type))
    return(edgeRenderInfo(Graph(object))$label)
  
  Afix <- index(object)$A ## Matrix with fixed parameters and ones where parameters are free
  Pfix <- index(object)$P ## Matrix with fixed covariance parameters and ones where param
  npar.mean <- index(object)$npar.mean
  Par <- object$coef
  if (npar.mean>0)
    Par <- Par[-(1:npar.mean),,drop=FALSE]
  Par <-
    switch(type,
           sd = paste(formatC(Par[,1,drop=FALSE]), " (", formatC(Par[,2,drop=FALSE]), ")", sep=""),
           est = formatC(Par[,1,drop=FALSE]),
           pval = formatC(Par[,4,drop=FALSE]),
           name = rownames(Par),
           none = ""
           )
  AP <- matrices(Model(object), Par)
  A <- AP$A; P <- AP$P
  P[exogenous(object),exogenous(object)] <- NA
  
  gr <- finalize(Model(object), diag=TRUE, cor=TRUE, addcolor=TRUE)
  Anz <- A; Anz[Afix==0] <- NA
  gr <- edgelabels(gr, lab=Anz)
  Pnz <- P; Pnz[Model(object)$cov==0] <- NA
  gr <- edgelabels(gr, lab=Pnz)
  Graph(object) <- gr
  return(object)
}

`edgelabels` <- function(object, ...) UseMethod("edgelabels")

`edgelabels<-` <- function(object,...,value) UseMethod("edgelabels<-")

`edgelabels<-.lvm` <- function(object,to,...,value) {
  edgelabels(object,to=to, lab=value,...)
}

`edgelabels<-.graphNEL` <- function(object,...,value) {
  edgelabels(object,lab=value,...)
}

`edgelabels.graphNEL` <- function(object, lab=NULL, to=NULL, from=NULL, cex=1.5, lwd=1, col="black", labcol="black",
                                  debug=FALSE,...) {
  if (is.null(lab)) {
    return(edgeRenderInfo(object)$label)
  }
  if (class(to)[1]=="formula") {
    yy <- getoutcome(to)
    from <- setdiff(all.vars(to),yy)
    to <- yy     
  }

  M <- as(object, Class="matrix")
  nodes <- nodes(object)

  if (is.null(edgeRenderInfo(object)$label))
    edgeRenderInfo(object)$label <- expression()
  
  if (!is.null(lab)) {
    if (!is.null(from) & !is.null(to)) {
      estr <- paste("\"",from,"~",to,"\"", sep="")
      estr2 <- paste(from,"~",to, sep="")
      Debug(estr,debug)
      if (!is.null(lab))
        edgeRenderInfo(object)$label[estr2] <- lab
      if (!is.null(cex))
        edgeRenderInfo(object)$cex[estr2] <- cex
      if (!is.null(col))
        edgeRenderInfo(object)$col[estr2] <- col
      if (!is.null(lwd))
        edgeRenderInfo(object)$lwd[estr2] <- lwd
      if (!is.null(labcol))
        edgeRenderInfo(object)$textCol[estr2] <- labcol
      st <- eval(parse(text=paste("expression(",edgeRenderInfo(object)$label,")",sep="")))
   ##   edgeRenderInfo(object)$label <- st
##      edgeRenderInfo(object)$label <- as.expression(edgeRenderInfo(object)$label)
      return(object)
    } 

    ## Used by "edgelabels.lvmfit"
    for (r in 1:nrow(M))
      for (s in 1:ncol(M)) {
        if (M[r,s]!=0 & !is.na(lab[r,s])) {
          estr <- paste("\"",nodes[r],"~",nodes[s],"\"", sep="")
          estr2 <- paste(nodes[r],"~",nodes[s], sep="")          
          Debug(estr, debug)
          st <- eval(parse(text=paste("expression(",lab[r,s],")",sep="")))
          edgeRenderInfo(object)$label[estr2] <- st
        }
      }
  }
  
  return(object)
}

`edgelabels.lvm` <-
function(object,lab=NULL,...) {
  if (is.null(lab)) {
    return(edgelabels(Graph(object),lab=lab,...))
  }
  Graph(object) <- edgelabels(Graph(object),lab=lab,...)
  return(object)
}

###}}} edgelabels

