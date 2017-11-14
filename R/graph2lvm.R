##' @export
`graph2lvm` <-
function(g, debug=FALSE, messages=0) {
  res <- lvm(graph::nodes(g), debug=debug, messages=messages)
  M <- t(as(g, Class="matrix"))
  for (i in seq_len(nrow(M))) {
    if (any(M[,i]==1)) {
      res <- regression(res, rownames(M)[M[,i]==1], rownames(M)[i], messages=messages)
    }
  }
  res
}

