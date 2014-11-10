##' @export
`measurement` <-
function(x, silent=TRUE, debug=FALSE) {
  M <- x$M
  latent.idx <- match(latent(x),vars(x))
  obs.idx <- match(manifest(x),vars(x))
  if (length(latent.idx)==0)
    return(NULL)
  
  measurementmodels <- c()
  for (i in 1:length(latent.idx)) {
    ii <- latent.idx[i]
    
    relation <- M[obs.idx,ii]==1
    byNodes <- names(relation)[relation]
    newnodes <- c(latent(x)[i],byNodes)
    g0 <- graph::subGraph(newnodes, Graph(x))
    lvm1 <- latent(graph2lvm(g0, debug=TRUE), latent(x)[i])
    g0fix<- x$fix[newnodes, newnodes]; lvm1$fix <- g0fix
    index(lvm1) <- reindex(lvm1)
    measurementmodels <- c(measurementmodels, list(lvm1))
  }

  measurementmodels    
}

