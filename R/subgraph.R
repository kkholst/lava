subgraph <- function(g,from,to,Tree=new("graphNEL",node=c(to,from),edgemode="directed"),...) {
  adjnodes <- adj(g,from)[[1]]
  newnodes <- !(adjnodes %in% nodes(Tree))
  if (length(adjnodes)==0)
    return(Tree)
  for (v in adjnodes) {  
    if (v==to) {
      Tree <- addEdge(from, v, Tree)
    }
    re1 <- acc(g,v)[[1]] ## Reachable nodes from v
    if ((to %in% names(re1)[re1>0])) {
      if (!(v %in% nodes(Tree)))
        Tree <- addNode(v,Tree)
      Tree <- addEdge(from, v, Tree)
      Tree <- path(g,v,to,Tree)
    }
  }
  return(Tree)
}

