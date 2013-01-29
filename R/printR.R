##' @title print vector/matrix in R format
##' @aliases tovec
##' @param x 
##' @param eol 
##' @param ... 
##' @author Klaus K. Holst
##' @export
printR <- function(x,eol="\n",...) {
  if (is.vector(x)) {
    row <- paste("c(",paste(x,collapse=","),")",sep="")
    cat(row,eol,sep="")
  }
  if (is.matrix(x)) {
    row <- "rbind("
    cat(row,eol)
    for (i in 1:nrow(x)) {
      printR(x[i,],"")
      cat(ifelse (i<nrow(x),",",""),eol,sep="")      
    }
    row <- ")"
    cat(row,eol)
  }
  invisible(x)
}
tovec <- function(x,...) paste("c(",paste(x,collapse=", "),")")
