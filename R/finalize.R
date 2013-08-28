##' @export
`finalize` <-
function(x,...) UseMethod("finalize")

##' @S3method finalize lvm
`finalize.lvm` <-
function(x, diag=FALSE, cor=FALSE, addcolor=TRUE, intercept=FALSE, plain=FALSE, cex, fontsize1=10, cols=c("lightblue","orange","yellowgreen"), unexpr=FALSE, addstyle=TRUE, ...) {

  g <- as(new("graphAM",adjMat=x$M,"directed"),"graphNEL")
  nodeRenderInfo(g)$fill <- NA
  nodeRenderInfo(g)$label <- NA
  nodeRenderInfo(g)$label[vars(x)] <- vars(x)
    
    nodeRenderInfo(g)$shape <- x$graphdef$shape
    Lab <- NULL
    for (i in seq_len(length(x$noderender))) {
      nn <- unlist(x$noderender[[i]])
      if (length(nn)>0) {
        R <- list(as.list(x$noderender[[i]])); names(R) <- names(x$noderender)[i]        
        if (names(x$noderender)[i]!="label")
          nodeRenderInfo(g) <- x$noderender[i]
        else Lab <- R[[1]]
      }
    }
    
    if (!is.null(Lab)) { ## Ugly hack to allow mathematical annotation
      nn <- names(nodeRenderInfo(g)$label)
      LL <- as.list(nodeRenderInfo(g)$label)
      LL[names(Lab)] <- Lab
      nodeRenderInfo(g) <- list(label=as.expression(LL))
      names(nodeRenderInfo(g)$label) <- nn
      ii <- which(names(nodeRenderInfo(g)$label)=="")
      if (length(ii)>0)
      nodeRenderInfo(g)$label <- nodeRenderInfo(g)$label[-ii]
    }
   
    edgeDataDefaults(g)$futureinfo <- x$edgerender$futureinfo    
    edgeRenderInfo(g)$lty <- x$graphdef$lty
    edgeRenderInfo(g)$lwd <- x$graphdef$lty
    edgeRenderInfo(g)$col <- x$graphdef$col
    edgeRenderInfo(g)$textCol <- x$graphdef$textCol
    edgeRenderInfo(g)$arrowhead <- x$graphdef$arrowhead
    edgeRenderInfo(g)$dir <- x$graphdef$dir
    edgeRenderInfo(g)$arrowtail <- "none"
    edgeRenderInfo(g)$cex <- x$graphdef$cex
    edgeRenderInfo(g)$label <- x$graphdef$label
    for (i in seq_len(length(x$edgerender))) {
      ee <- x$edgerender[[i]]
      if (length(ee)>0 && names(x$edgerender)[i]!="futureinfo") {
        edgeRenderInfo(g)[names(x$edgerender)[i]][names(ee)] <- ee
      }
    }
  
  opt <- options(warn=-1)
  var <- rownames(covariance(x)$rel)


   if (unexpr) {
    mylab <- as.character(edgeRenderInfo(g)$label); names(mylab) <- names(edgeRenderInfo(g)$label)
    g@renderInfo@edges$label <- as.list(mylab)
  } 
 
  
  if (intercept) {
  ##  mu <- intfix(x)
  ##  nNA <- sum(is.na(mu))
 ##   if (nNA>0)
##      mu[is.na(mu)] <- paste("m",1:nNA)
##    mu <- unlist(mu)
##    x <- addNode(mu,x)
##    for (i in 1:length(mu)) {
  ##    print(mu[i])
##      x <- addEdge(var[i], var[i], x)
##    }
##    x <- addattr(x,attr="shape",var=mu,val="none")      
  }

  allEdges <- edgeNames(g)
  regEdges <- c()
  feedback <- c()
  A <- index(x)$A
  if (index(x)$npar.reg>0)
  for (i in 1:(nrow(A)-1))
      for (j in (i+1):(ncol(A))) {
      if(A[i,j]==1 & A[j,i]==1) feedback <- c(feedback,
                          paste(var[i],"~",var[j], sep=""),
                          paste(var[j],"~",var[i], sep=""))
      if (A[j,i]==0 & x$M[j,i]!=0) {
          g <- removeEdge(var[j],var[i],g)
      }
      if (A[i,j]==1) regEdges <- c(regEdges,paste(var[i],"~",var[j], sep=""))
      if (A[j,i]==1) regEdges <- c(regEdges,paste(var[j],"~",var[i], sep=""))
    }

  
  varEdges <- corEdges <- c()
  delta <- ifelse(diag,0,1)  
  if (cor | diag) {
    for (r in 1:(nrow(covariance(x)$rel)-delta) ) {
      for (s in (r+delta):ncol(covariance(x)$rel) ) {
        if (cor | r==s)
          if (covariance(x)$rel[r,s]==1 & (!any(c(var[r],var[s])%in%exogenous(x)))) {
            newedges <- c()
            if (A[r,s]!=1) {
              g <- addEdge(var[r],var[s], g)
              newedges <- paste(var[r],"~",var[s], sep="")
            } else {
              if (A[s,r]!=1) {
                g <- addEdge(var[s],var[r], g)
                newedges <- c(newedges,paste(var[s],"~",var[r], sep=""))
              }
            }
##            g <- addEdge(var[s],var[r], g)
##            newedges <- c(paste(var[r],"~",var[s], sep=""),
##                          paste(var[s],"~",var[r], sep=""))
            if (r==s)
              varEdges <- c(varEdges,
                            newedges
                            )
            if (r!=s)
            corEdges <- c(corEdges,newedges)
          }
      }
    }
  }

  if (length(x$edgerender$futureinfo)>0) {
    estr <- names(x$edgerender$futureinfo$label)
    estr <- estr[which(unlist(lapply(estr,nchar))>0)]
    revestr <- sapply(estr, function(y) paste(rev(unlist(strsplit(y,"~"))),collapse="~"))
    revidx <- which(revestr%in%edgeNames(g))
    count <- 0
    for (i in estr) {      
      count <- count+1
      for (f in names(x$edgerender$futureinfo)) {
        if (count%in%revidx) {
          g@renderInfo@edges[[f]][[revestr[count]]] <- x$edgerender$futureinfo[[f]][[i]]
        } else {
          g@renderInfo@edges[[f]][[i]] <- x$edgerender$futureinfo[[f]][[i]]
        }
      }
    }
  }
  allEdges <- unique(c(regEdges,corEdges,varEdges))
  corEdges <- setdiff(corEdges,regEdges)

  for (e in allEdges) {
    dir <- "forward"; lty <- 1; arrowtail <- "none"
    if (e %in% feedback) {
      dir <- "none"; lty <- 1; arrowtail <- "open"
    }
    if (e %in% varEdges) {
      dir <- "none"; lty <- 2; arrowtail <- "none"
    }
    if (e %in% corEdges) {
      dir <- "none"; lty <- 2; arrowtail <- "open"
    }
    arrowhead <- "open"
##    estr <- paste("\"",e,"\"",sep="")
    estr <- e
    for (f in c("col","cex","textCol","lwd","lty")) {
      if (!(estr%in%names(edgeRenderInfo(g)[[f]]))
          || is.na(edgeRenderInfo(g)[[f]][[estr]]))
        g <- addattr(g,f,var=estr,
                     val=x$graphdef[[f]],
                     fun="edgeRenderInfo")      
    }


    if (addstyle) {
      g <- addattr(g,"lty",var=estr,val=lty,fun="edgeRenderInfo")
      g <- addattr(g,"direction",var=estr,val=dir,fun="edgeRenderInfo")
      g <- addattr(g,"dir",var=estr,val=dir,fun="edgeRenderInfo")
      g <- addattr(g,"arrowhead",var=estr,val=arrowhead,fun="edgeRenderInfo")
      g <- addattr(g,"arrowtail",var=estr,val=arrowtail,fun="edgeRenderInfo")
      g <- addattr(g,attr="fontsize",var=estr,val=fontsize1,fun="edgeRenderInfo")
    }
    if (is.null(edgeRenderInfo(g)$label))
      edgeRenderInfo(g)$label <- expression()

    if (!missing(cex))
      if (!is.null(cex))
        nodeRenderInfo(g)$cex <- cex
  }
  if (plain) {
    g <- addattr(g,attr="shape",var=vars(x),val="none")      
  } else {    
    if (addcolor) {
      if (is.null(x$noderender$fill)) notcolored <- vars(x)
      else notcolored <- vars(x)[is.na(x$noderender$fill)]
      nodecolor(g, intersect(notcolored,exogenous(x))) <- cols[1]
      nodecolor(g, intersect(notcolored,endogenous(x))) <- cols[2]
      nodecolor(g, intersect(notcolored,latent(x))) <- cols[3]
      ##        nodecolor(x, intersect(notcolored,survival(x))) <- cols[4]        
      myhooks <- gethook("color.hooks")
      count <- 3
      for (f in myhooks) {
        count <- count+1
        res <- do.call(f, list(x=x,subset=notcolored))
        if (length(cols)>=count) {
          nodecolor(g,res$vars) <- cols[count]
        } else {
          nodecolor(g, res$vars) <- res$col
        }
      }       
    }
  }
  options(opt)
  attributes(g)$feedback <- (length(feedback)>0)
  return(g)
}


