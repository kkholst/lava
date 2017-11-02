
impute0 <- function(object,rows,idx,na.action=na.omit,...) {
    if (missing(rows) && missing(idx)) {
        df <- na.action(object,...) 
        rows <- attr(df,"na.action")  
    }
    if (!missing(idx)) {
        obs1 <- setdiff(seq(length(object)),idx)[1]
    } else {
        obs1 <- setdiff(seq(NROW(object)),rows)[1]
    }
    fobs <- object[obs1]
    if (is.logical(fobs)) val <- FALSE
    else if (is.character(fobs)) val <- fobs
    else if (is.factor(fobs)) val <- levels(fobs)[1]
    else val <- 0
    if (!missing(idx)) {
        object[idx] <- val
        return(object)
    }
    if (is.matrix(object)) {
        object[rows,] <- val
    } else {
        object[rows] <- val
    }
    return(object)
}

##' @export 
na.pass0 <- function(object,all=TRUE,na.action=na.omit, ...) {
    ## Fill in "zeros" in the design matrix where we have missing data    
    df <- na.action(object,...)
    idx <- attr(df,"na.action")
    if (is.matrix(object) || is.vector(object)) {
        object <- impute0(object,rows=idx)
    } else {
        for (i in seq_len(NCOL(object))) {
            object[[i]] <- impute0(object[[i]],rows=idx)
        }
    }
    structure(object,na.action=structure(idx,class="pass0"))
}
