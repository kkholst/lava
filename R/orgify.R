##' @title Convert object to ascii suitable for org-mode
##' @param x 
##' @param ncol 
##' @param include.rownames 
##' @param include.colnames 
##' @param header 
##' @param frame 
##' @param rownames 
##' @param type 
##' @param ... 
##' @author Klaus K. Holst
##' @export
orgify <- function(x,ncol,include.rownames=TRUE,include.colnames=TRUE,header=TRUE, frame="topbot",rownames,type="org",...) {
  if (!missing(ncol)) {
    y <- formatC(as.vector(x))  
    n0 <- length(y)%%ncol
    if (n0 > 0) 
        y <- c(y, rep("", ncol - n0))
    res <- ascii(matrix(y, ncol = ncol, byrow = TRUE),type=type,...)
    return(res)
  }
  if (is.vector(x)) {
    if (is.null(names(x))) {
      include.colnames <- FALSE
      header <- FALSE
    }
    x <- rbind(x)
    if (!missing(rownames)) {
      rownames(x) <- rownames[1]
    } else {
      include.rownames <- FALSE
    }      
  }
  res <- ascii(x,include.rownames=include.rownames,include.colnames=include.colnames,header=header,frame=frame,type=type,...)
  return(res)
}
