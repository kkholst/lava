##' @export
`addattr` <-
function(x,...) UseMethod("addattr")

##' @export
`addattr.lvm` <-
function(x, attr, var=NULL, val=TRUE, fun="graph::nodeRenderInfo",debug=FALSE,...) {
  if (!is.null(var)) {
    Graph(x) <- addattr(Graph(x), attr=attr, var=var, val=val, fun=fun, debug=debug)
    return(x)
  } else {
    addattr(Graph(x), attr=attr, var=var, val=val, fun=fun)
  }
}

##' @export
`addattr.graphNEL` <-
function(x, attr, var=NULL, val=TRUE,fun="graph::nodeRenderInfo",debug=FALSE,...) {
  if (is.null(var)) {
    f <- eval(call(fun,x))
    if (is.null(val))
      attrvar <- names(f[[attr]])
    else
      attrvar <- names(f[[attr]])[which(val==f[[attr]])]
    return(attrvar)
  }
  if (is.character(val))
    myexpr <- paste0("list(",attr,"=c(", paste0("\"",var,"\"=\"",val,"\"" , collapse=", "), "))")
  else
    myexpr <- paste0("list(",attr,"=c(", paste0("\"",var,"\"=",val, collapse=", "), "))")
  Debug(list("str=",myexpr),debug)
  eval(parse(text=paste0(fun,"(x) <- ",myexpr)))
  return(x)
}
