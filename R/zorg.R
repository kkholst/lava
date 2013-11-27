##' Convert object to ascii suitable for org-mode
##' 
##' @title Convert object to ascii suitable for org-mode
##' @param x R object
##' @param ... 
##' @param ncol If \code{x} is a vector and \code{ncol} is given as argument, the resulting output will be a \code{matrix} with \code{ncol} columns
##' @param include.rownames If \code{FALSE} row names are removed
##' @param include.colnames If \code{FALSE} column names are removed
##' @param header If TRUE the header is included
##' @param frame Frame argument (see \code{ascii})
##' @param rownames Optional vector of row names
##' @param colnames Optional vector of column names
##' @param type Type argument (see \code{ascii})
##' @param tab Tabulate?
##' @param margins Add margins to table?
##' @param print print or return result 
##' @param html HTML prefix (added to ATTR_HTML)
##' @param latex LaTeX prefix (added to ATTR_LaTeX)
##' @param \dots additional arguments to lower level functions
##' @author Klaus K. Holst
org <- function(x,...,ncol,include.rownames=TRUE,include.colnames=TRUE,header=TRUE, frame="topbot",rownames=NULL,colnames=NULL,type="org",tab=FALSE,margins=TRUE,print=TRUE,html,latex) {  
    if (!suppressPackageStartupMessages(require(ascii))) stop("ascii package required")
    if (tab) {
        if (!inherits(x,"table")) {
            x <- table(x)
        }
        if (margins) x <- addmargins(x)
    }
    if (!missing(ncol)) {
        y <- formatC(as.vector(x))  
        n0 <- length(y)%%ncol
        if (n0 > 0) 
            y <- c(y, rep("", ncol - n0))
        x <- matrix(y, ncol = ncol, byrow = TRUE)
    }
    if (is.vector(x)) {
        if (is.null(names(x))) {
            include.colnames <- FALSE
            header <- FALSE
        }
        x <- rbind(x)
        if (!is.null(rownames)) {
            rownames(x) <- rownames[1]
        } else {
            include.rownames <- FALSE
        }      
    }
    x <- ascii(x,include.rownames=include.rownames,include.colnames=include.colnames,header=header,frame=frame,type=type,rownames=rownames,colnames=colnames,...)
    if (print) {
        op <- options(asciiType=type)
        if (!missing(html)) 
            cat("#+ATTR_HTML: ",html,"\n",sep="")            
        if (!missing(latex)) 
            cat("#+ATTR_LaTeX: ",latex,"\n",sep="")            
        suppressWarnings(ascii::print(x,...))
        options(op)
    }
    invisible(x)
}

