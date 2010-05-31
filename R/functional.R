"functional<-" <- function(x,...,value) UseMethod("functional<-")
"functional<-.lvm" <- function(x,to,from,...,value) {
  if (class(to)[1]=="formula") {
    yy <- decomp.specials(getoutcome(to))
    ##xx <- attributes(terms(to))$term.labels
    myvars <- all.vars(to)
    xx <- setdiff(myvars,yy)
    if (length(yy)*length(xx)>length(value) & length(value)!=1) stop("Wrong number of values")
    count <- 0
    for (y in yy) {
      count <- count+1
      for (i in 1:length(xx)) {
        suppressWarnings(x <- regression(x, to=y,from=xx[i],silent=TRUE))
        count <- count+1
          if (length(value)==1) {
            functional(x, to=y, from=xx[i],...) <- value
          } else
          functional(x, to=y, from=xx[i],...) <- value[[count]]        
        }
    }
    return(x)
  }
  
  if (missing(from) | missing(to))
    return(x)

  edges <- paste(from,to,sep="~")
  edges. <- paste("\"", edges, "\"", sep="")
  mytext <- paste("c(", paste(paste(edges.,"=",expression(value),sep=""),collapse=","),")")
  edgeRenderInfo(Graph(x))$"functional" <- eval(parse(text=mytext))
  return(x)
}

"functional" <- function(x,...) UseMethod("functional")
functional.lvm <- function(x,to,from,...) {
  if (missing(from))
    return(edgeRenderInfo(Graph(x))$functional)
  
  edges <- paste(from,to,sep="~")
  edgeRenderInfo(Graph(x))$functional[edges]
}

