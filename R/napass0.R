
impute0 <- function(object,rows,idx,na.action=na.omit,value,...) {
  if (missing(rows) && missing(idx)) {
    df <- na.action(object,...)
    rows <- attr(df,"na.action")
  }
  if (!missing(idx)) {
    obs1 <- setdiff(seq(length(object)),idx)[1]
  } else {
    obs1 <- setdiff(seq(NROW(object)),rows)[1]
  }
  if (missing(value)) {
    fobs <- object[obs1]
    if (is.logical(fobs)) value <- FALSE
    else if (is.character(fobs)) value <- fobs
    else if (is.factor(fobs)) value <- levels(fobs)[1]
    else value <- 0
  }
  if (!missing(idx)) {
    object[idx] <- value
    return(object)
  }
  if (is.matrix(object)) {
    object[rows,] <- value
  } else {
    object[rows] <- value
  }
  return(object)
}

##' @title Handle Missing Values in Objects
##' @description Returns the object with missing data replaced by zeros. This is
##'   sometimes useful for example when working with inverse probability
##'   weighting of the complete-case data.
##' @seealso [na.pass()], [na.omit()], [na.fail()]
##' @export
##' @param object data.frame, vector
##' @param na.action function used to identify missing values
##' @param row.wise when true entire row is set to zero whenever a missing data
##'   is observed in that row
##' @param value value to replace NAs with (default 0)
##' @param ... additional arguments to lower level functions
##' @examples
##' d <- data.frame(y=c(1,1,NA,2,NA,2), r=c(1,1,0,1,1,1))
##' na.pass0(d)
##' glm(y ~ 1, weights=d$r, data=d, na.action=na.pass0)
na.pass0 <- function(object,
                     na.action = na.omit,
                     row.wise = FALSE,
                     value = 0, ...) {
  if (!row.wise && NCOL(object)>1L) {
    nas <- is.na(object)
    if (!any(nas)) return(object)
    idx <- which(nas)
    ii <- arrayInd(idx, dim(object))
    if (length(ii) > 0) {
      object[ii] <- value
    }
    return(structure(object,na.action=structure(idx,class="pass0")))
  }
  df <- na.action(object,...)
  idx <- attr(df, "na.action")
  if (is.null(idx)) return(object)
  if (is.matrix(object) || is.vector(object)) {
    object <- impute0(object,idx=idx,value=value,...)
  } else {
    for (i in seq_len(NCOL(object))) {
      object[[i]] <- impute0(object[[i]],idx=idx,value=value,...)
    }
  }
  return(structure(object,na.action=structure(idx,class="pass0")))
}
