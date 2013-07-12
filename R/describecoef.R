##' @export
describecoef <- function(x,par,from,to,mean=TRUE) {  
  p <- coef(x, mean=mean)  
  if (!missing(from)) {
    st1 <- paste(to,lava.options()$symbol[1],from,sep="")
    st2 <- paste(to,lava.options()$symbol[2],from,sep="")
    st3 <- paste(from,lava.options()$symbol[2],to,sep="")
    pos <- na.omit(match(unique(c(st1,st2,st3)),p))
    attributes(pos) <- NULL
    return(pos)
  }
  res <- strsplit(p,lava.options()$symbol[2])
  var.idx <- which(unlist(lapply(res,length))>1) ## Variance parameters
  rest.idx <- setdiff(1:length(p),var.idx)
  res[rest.idx] <- strsplit(p[rest.idx],lava.options()$symbol[1])
  mean.idx <- which(unlist(lapply(res,length))==1) ## Mean parameters
  reg.idx <- setdiff(rest.idx,mean.idx)
  names(res)[mean.idx] <- paste("m",1:length(mean.idx),sep="")
  for (i in var.idx)
    attr(res[[i]],"type") <- "cov"
  for (i in mean.idx)
    attr(res[[i]],"type") <- "mean"
  for (i in reg.idx)
    attr(res[[i]],"type") <- "reg"
  if (missing(par))
    return(res)
  return(res[par])
}

