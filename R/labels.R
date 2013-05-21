###{{{ labels

##' Define labels of graph
##' 
##' Alters labels of nodes and edges in the graph of a latent variable model
##' 
##' 
##' @aliases labels<- labels labels<-.default labels.lvm labels.lvmfit
##' labels.graphNEL edgelabels edgelabels<- edgelabels<-.lvm nodecolor
##' nodecolor<- nodecolor<-.default
##' @author Klaus K. Holst
##' @export
##' @keywords graphs aplot
##' @examples
##' m <- lvm(c(y,v)~x+z)
##' regression(m) <- c(v,x)~z
##' labels(m) <- c(y=expression(psi), z=expression(zeta))
##' nodecolor(m,~y+z+x,border=c("white","white","black"),
##'           labcol="white", lwd=c(1,1,5)) <-  c("orange","indianred","lightgreen")
##' edgelabels(m,y~z+x, cex=c(2,3), col=c("orange","black"),labcol="darkblue",
##'            lwd=c(3,1)) <- expression(phi,rho)
##' edgelabels(m,c(v,x)~z, labcol="red", cex=2) <- 2
##' \donttest{plot(m)}
##' @param object \code{lvm}-object.
##' @param value node label/edge label/color
##' @param to Formula specifying outcomes and predictors defining relevant
##' edges.
##' @param \dots Additional arguments (\code{lwd}, \code{cex}, \code{col},
##' \code{labcol}), \code{border}.
##' @param var Formula or character vector specifying the nodes/variables to
##' alter.
##' @param border Colors of borders
##' @param labcol Text label colors
##' @param shape Shape of node
##' @param lwd Line width of border
##' @usage
##' \method{labels}{default}(object, ...) <- value
##' \method{edgelabels}{lvm}(object, to, ...) <- value
##' \method{nodecolor}{default}(object, var=vars(object),
##' border, labcol, shape, lwd, ...) <- value
`labels<-` <- function(object,...,value) UseMethod("labels<-")

##' @S3method labels<- default
`labels<-.default` <- function(object,...,value) {
  labels(object,value)
}

##' @S3method labels graphNEL
labels.graphNEL <- function(object,lab=NULL,...) {
  if (is.null(lab))    
    return(nodeRenderInfo(object)$label)  
  nodeRenderInfo(object) <- list(label=lab)
  names(nodeRenderInfo(object)$label) <- nodes(object);
  return(object)
}

##' @S3method labels lvmfit
labels.lvmfit <- function(object,lab=NULL,...) {
  if (is.null(lab)) return(object$noderender$label)
  object$noderender$label <- lab
  return(object)
}

##' @S3method labels lvm
`labels.lvm` <- function(object,lab=NULL,...) {
  if (is.null(lab))
    return(object$noderender$label)
  if (is.null(object$noderender$label))
    object$noderender$label <- lab
  else
    object$noderender$label[names(lab)] <- lab
  return(object)
}
###}}} labels

###{{{ edgelabels

##' @S3method edgelabels<- lvmfit
"edgelabels<-.lvmfit" <- function(object,to,from,est=TRUE,edges=NULL,cex=1,...,value) {
  if (is.null(edges))  {
    if (class(to)[1]=="formula") {
      yy <- decomp.specials(getoutcome(to))
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

##' @S3method edgelabels lvmfit
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
  
  gr <- finalize(Model(object), ...)
  Anz <- A; Anz[Afix==0] <- NA
  gr <- edgelabels(gr, lab=Anz)
  Pnz <- P; Pnz[Model(object)$cov==0] <- NA
  gr <- edgelabels(gr, lab=Pnz)
  Graph(object) <- gr
  return(object)
}

##' @export
`edgelabels` <- function(object, ...) UseMethod("edgelabels")

##' @export
`edgelabels<-` <- function(object,...,value) UseMethod("edgelabels<-")

##' @S3method edgelabels<- lvm
`edgelabels<-.lvm` <- function(object,to,...,value) {
  edgelabels(object,to=to, lab=value,...)
}

##' @S3method edgelabels<- graphNEL
`edgelabels<-.graphNEL` <- function(object,...,value) {
  edgelabels(object,lab=value,...)
}

##' @S3method edgelabels graphNEL
`edgelabels.graphNEL` <- function(object, lab=NULL, to=NULL, from=NULL, cex=1.5, lwd=1, lty=1, col="black", labcol="black",
                                  expr=TRUE,
                                  debug=FALSE,...) {
  if (is.null(lab)) {
    return(edgeRenderInfo(object)$label)
  }
  if (class(to)[1]=="formula") {
    yy <- decomp.specials(getoutcome(to))
    from <- all.vars(to[[3]])##setdiff(all.vars(to),yy)
    if (length(from)==0) from <- yy
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
      if (length(lab)!=length(estr2)) lab <- rep(lab,length(estr2))
      if (length(col)!=length(estr2)) col <- rep(col,length(estr2))
      if (length(cex)!=length(estr2)) cex <- rep(cex,length(estr2))
      if (length(lwd)!=length(estr2)) lwd <- rep(lwd,length(estr2))
      if (length(lty)!=length(estr2)) lty <- rep(lty,length(estr2))
      if (length(labcol)!=length(estr2)) labcol <- rep(labcol,length(estr2))  

      curedges <- names(edgeRenderInfo(object)$label)
       Debug(estr,debug)
      
      estr2.idx <- which(estr2%in%curedges)
      newstr.idx <- setdiff(1:length(estr2),estr2.idx)
      newstr <- estr2[newstr.idx]
      estr2 <- estr2[estr2.idx]
      if (length(estr2)>0) {
        if (!is.null(lab))
          edgeRenderInfo(object)$label[estr2] <- lab[estr2.idx]
        if (!is.null(cex))
        edgeRenderInfo(object)$cex[estr2] <- cex[estr2.idx]
        if (!is.null(col))
          edgeRenderInfo(object)$col[estr2] <- col[estr2.idx]
        if (!is.null(lwd))
          edgeRenderInfo(object)$lwd[estr2] <- lwd[estr2.idx]
        if (!is.null(lty))
          edgeRenderInfo(object)$lty[estr2] <- lty[estr2.idx]
        if (!is.null(labcol))
          edgeRenderInfo(object)$textCol[estr2] <- labcol[estr2.idx]
      }
      if (length(newstr)>0) {        
        
        if (!is.null(lab))
          edgeDataDefaults(object)$futureinfo$label[newstr] <-
            lab[newstr.idx]
        if (!is.null(cex))
          edgeDataDefaults(object)$futureinfo$cex[newstr] <-
            cex[newstr.idx]
        if (!is.null(col))
          edgeDataDefaults(object)$futureinfo$col[newstr] <-
            col[newstr.idx]
        if (!is.null(lwd))
          edgeDataDefaults(object)$futureinfo$lwd[newstr] <-
            lwd[newstr.idx]
        if (!is.null(lty))
          edgeDataDefaults(object)$futureinfo$lty[newstr] <-
            lty[newstr.idx]
        if (!is.null(labcol))        
          edgeDataDefaults(object)$futureinfo$textCol[newstr] <-
            labcol[newstr.idx]
      }
      return(object)
    } 

    ## Used by "edgelabels.lvmfit"
    for (r in 1:nrow(M))
      for (s in 1:ncol(M)) {
        if (M[r,s]!=0 & !is.na(lab[r,s])) {
          estr <- paste("\"",nodes[r],"~",nodes[s],"\"", sep="")
          estr2 <- paste(nodes[r],"~",nodes[s], sep="")          
          Debug(estr, debug)
          if (expr)
            st <- eval(parse(text=paste("expression(",lab[r,s],")",sep="")))
          else
            st <- lab[r,s]
          edgeRenderInfo(object)$label[estr2] <- st
        }
      }
  }
  
  return(object)
}



##' @S3method edgelabels lvm
`edgelabels.lvm` <- function(object, lab=NULL, to=NULL, from=NULL,
                             cex=1.5, lwd=1, lty=1, col="black", labcol="black",
                             expr=TRUE, debug=FALSE,...) {
  if (is.null(lab)) {
    return(object$edgerender$label)
  }
  if (class(to)[1]=="formula") {
    yy <- decomp.specials(getoutcome(to))
    from <- all.vars(to[[3]])##setdiff(all.vars(to),yy)
    if (length(from)==0) from <- yy
    to <- yy     
  }

  M <- object$M
  nodes <- colnames(M)

  if (is.null(object$edgerender$label))
    object$edgerender$label <- expression()

  
  if (!is.null(lab)) {
    if (!is.null(from) & !is.null(to)) {      
      estr <- paste("\"",from,"~",to,"\"", sep="")
      estr2 <- paste(from,"~",to, sep="")
      if (length(lab)!=length(estr2)) lab <- rep(lab,length(estr2))
      if (length(col)!=length(estr2)) col <- rep(col,length(estr2))
      if (length(cex)!=length(estr2)) cex <- rep(cex,length(estr2))
      if (length(lwd)!=length(estr2)) lwd <- rep(lwd,length(estr2))
      if (length(lty)!=length(estr2)) lty <- rep(lty,length(estr2))
      if (length(labcol)!=length(estr2)) labcol <- rep(labcol,length(estr2))  

      curedges <- names(object$edgerender$label)
       Debug(estr,debug)
      
      estr2.idx <- which(estr2%in%curedges)
      newstr.idx <- setdiff(1:length(estr2),estr2.idx)
      newstr <- estr2[newstr.idx]
      estr2 <- estr2[estr2.idx]
      if (length(estr2)>0) {
        if (!is.null(lab))
          object$edgerenderlabel[estr2] <- lab[estr2.idx]
        if (!is.null(cex))
          object$edgerender$cex[estr2] <- cex[estr2.idx]
        if (!is.null(col))
          object$edgerender$col[estr2] <- col[estr2.idx]
        if (!is.null(lwd))
          object$edgerender$lwd[estr2] <- lwd[estr2.idx]
        if (!is.null(lty))
          object$edgerender$lty[estr2] <- lty[estr2.idx]
        if (!is.null(labcol))
          object$edgerender$textCol[estr2] <- labcol[estr2.idx]
      }
      if (length(newstr)>0) {        
        
        if (!is.null(lab))
          object$edgerender$futureinfo$label[newstr] <-
            lab[newstr.idx]
        if (!is.null(cex))
          object$edgerender$futureinfo$cex[newstr] <-
            cex[newstr.idx]
        if (!is.null(col))
          object$edgerender$futureinfo$col[newstr] <-
            col[newstr.idx]
        if (!is.null(lwd))
          object$edgerender$futureinfo$lwd[newstr] <-
            lwd[newstr.idx]
        if (!is.null(lty))
          object$edgerender$futureinfo$lty[newstr] <-
            lty[newstr.idx]
        if (!is.null(labcol))        
          object$edgerender$futureinfo$textCol[newstr] <-
            labcol[newstr.idx]
      }
      return(object)
    } 

    ## Used by "edgelabels.lvmfit"
    for (r in 1:nrow(M))
      for (s in 1:ncol(M)) {
        if (M[r,s]!=0 & !is.na(lab[r,s])) {
          estr <- paste("\"",nodes[r],"~",nodes[s],"\"", sep="")
          estr2 <- paste(nodes[r],"~",nodes[s], sep="")          
          Debug(estr, debug)
          if (expr)
            st <- eval(parse(text=paste("expression(",lab[r,s],")",sep="")))
          else
            st <- lab[r,s]
          object$edgerender$label[estr2] <- st
        }
      }
  }  
  return(object)
}

###}}} edgelabels

