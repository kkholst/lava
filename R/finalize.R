`finalize` <-
function(model,...) UseMethod("finalize")

`finalize.lvm` <-
function(model, diag=FALSE, cor=FALSE, addcolor=TRUE, intercept=FALSE, plain=FALSE, cex, fontsize1=10, cols=c("lightblue","orange","yellowgreen","indianred1"), ...) {
  opt <- options(warn=-1)
  x <- model$graph
  
  var <- rownames(covariance(model)$rel)
  edges <- c()
  delta <- ifelse(diag,0,1)
  
  if (cor | diag) {
    for (r in 1:(nrow(covariance(model)$rel)-delta) ) {
      for (s in (r+delta):ncol(covariance(model)$rel) ) {
        if (cor | r==s)
          if (covariance(model)$rel[r,s]==1 & (!any(c(var[r],var[s])%in%exogenous(model)))) {
            x <- addEdge(var[r],var[s], x)
            x <- addEdge(var[s],var[r], x)
            edges <- c(edges,
                       paste(var[r],"~",var[s], sep=""),
                       paste(var[s],"~",var[r], sep=""))
          }
      }
    }
  }
  if (intercept) {
  ##  mu <- intfix(model)
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

  
  recursive <- c()

  A <- index(model)$A
  if (index(model)$npar.reg>0)
  for (i in 1:(nrow(A)-1))
    for (j in (i+1):(ncol(A)))
      if(A[i,j]==1 & A[j,i]==1) recursive <- c(recursive,
                        paste(var[i],"~",var[j], sep=""),
                        paste(var[j],"~",var[i], sep=""))
 
  for (e in edgeNames(x)) {
    dir <- "forward"; lty <- 1; arrowhead <- "none"
    if (e %in% edges) {
      dir <- "none"; lty <- 2; arrowhead <- "open"
    }
    if (e %in% recursive) {
      dir <- "none"; lty <- 1; arrowhead <- "open"
    }    
##    estr <- paste("\"",e,"\"",sep="")
    estr <- e
    x <- addattr(x,"dir",var=estr,val=dir,fun="edgeRenderInfo")
    x <- addattr(x,"direction",var=estr,val=dir,fun="edgeRenderInfo")
    x <- addattr(x,"arrowhead",var=estr,val="open",fun="edgeRenderInfo")
    x <- addattr(x,"arrowtail",var=estr,val=arrowhead,fun="edgeRenderInfo")
    x <- addattr(x,attr="lty",var=estr,val=lty,fun="edgeRenderInfo")
    x <- addattr(x,attr="fontsize",var=estr,val=fontsize1,fun="edgeRenderInfo")
    if (is.null(edgeRenderInfo(x)$label))
      edgeRenderInfo(x)$label <- expression()

    if (!missing(cex))
      if (!is.null(cex))
        nodeRenderInfo(x)$cex <- cex

    if (plain) {
      x <- addattr(x,attr="shape",var=vars(model),val="none")      
    } else {    
      if (addcolor) {
        if (is.null(nodeRenderInfo(Graph(model))$fill)) notcolored <- vars(model)
        else notcolored <- vars(model)[is.na(nodeRenderInfo(Graph(model))$fill)]
        nodecolor(x, intersect(notcolored,exogenous(model))) <- cols[1]
        nodecolor(x, intersect(notcolored,endogenous(model))) <- cols[2]
        nodecolor(x, intersect(notcolored,latent(model))) <- cols[3]
        nodecolor(x, intersect(notcolored,survival(model))) <- cols[4]
        nodecolor(x, intersect(notcolored,binary(model))) <- cols[4]
      }
    }
    
    
  }
  
  options(opt)
  return(x)
}


