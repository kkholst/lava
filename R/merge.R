##' @export
`%++%.lvm` <- function(x, y) merge(x, y)

##' @export
"+.lvm" <- function(x, ...) {
  merge(x, ...)
}

##' @export
merge.lvm <- function(x, y, ...) {
  objects <- list(x, y, ...)
  if (length(objects)<2) return(x)
  m <- objects[[1]]
  for (i in seq(2, length(objects))) {
    m2 <- objects[[i]]
    if (length(latent(m2))>0)
      latent(m) <- latent(m2)
    if (length(m2$constrain)>0)
      m$constrain <- c(m$constrain,m2$constrain)
    M <- (index(m2)$M)
    P <- (index(m2)$P)
    nn <- vars(m2)
    for (j in seq_len(nrow(M))) {
      if (any(idx <- M[j,]!=0)) {
        val <- as.list(rep(NA,sum(idx==TRUE)))
        if (any(idx. <- !is.na(m2$par[j,idx])))
          val[idx.] <- m2$par[j,idx][idx.]
        if (any(idx. <- !is.na(m2$fix[j,idx])))
          val[idx.] <- m2$fix[j,idx][idx.]
        regression(m,to=nn[idx],from=nn[j],messages=0) <- val
      }

      P0 <- P[j,]; P0[seq_len(j-1)] <- 0
      idx <- P[j,]!=0 | m2$covfix[j,]==0
      idx[is.na(idx)] <- FALSE
      if (any(idx)) {
          val <- as.list(rep(NA,sum(idx==TRUE)))
          if (any(idx. <- !is.na(m2$covpar[j,idx])))
            val[idx.] <- m2$covpar[j,idx][idx.]
          if (any(idx. <- !is.na(m2$covfix[j,idx])))
            val[idx.] <- m2$covfix[j,idx][idx.]
          covariance(m,nn[idx],nn[j],messages=0) <- val
      }
    }
    intercept(m,nn) <- intercept(m2)
    m2x <- exogenous(m2)
    if (length(m2x)>0)
      exogenous(m) <- c(exogenous(m),m2x)
  }
  index(m) <- reindex(m)
  return(m)
}

##' @export
merge.list <- function(x, ...) {
  do.call(merge, c(x, list(...)))
}

##' @export
merge.lm <- function(x,y,...) {
  args <- c(list(x, y), list(...))
  nn <- names(formals(merge.estimate)[-seq(3)])
  idx <- na.omit(match(nn, names(args)))
  models <- args
  models[idx] <- NULL
  mm <- lapply(models, function(x) tryCatch(estimate(x),error=function(e) NULL))
  names(mm)[1:2] <- c("x", "y")
  ii <- which(unlist(lapply(mm, is.null)))
  if (length(ii) > 0) mm[ii] <- NULL
  do.call(merge,c(mm,args[idx]))
}

##' @export
merge.glm <- merge.lm

##' @export
merge.lvmfit <- merge.lm

##' @export
merge.multinomial <- function(x,...) {
    merge.estimate(x,...)
}

