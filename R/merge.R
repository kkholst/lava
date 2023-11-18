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
"+.estimate" <- function(x, ...) {
  merge(x, ...)
}

##' @export
"-.estimate" <- function(x,...) {
  res <- merge(x, ...)
  estimate(res, pairwise.diff(length(coef(res))))
}

##' @export
merge.estimate <- function(x,y,...,
                           id,
                           paired=FALSE,
                           labels=NULL,
                           keep=NULL,
                           subset=NULL,
                           regex=FALSE,
                           ignore.case=FALSE) {
    objects <- list(x, estimate(y), ...)
    if (length(nai <- names(objects)=="NA")>0)
    names(objects)[which(nai)] <- ""
    if (!missing(subset)) {
      if (regex) {

      }
      coefs <- unlist(lapply(objects, function(x) coef(x)[subset]))
    } else {
      coefs <- unlist(lapply(objects,coef))
    }
    if (!is.null(labels)) {
      names(coefs) <- labels
    } else {
      names(coefs) <- make.unique(names(coefs))
    }
    if (regex) {
      if (!is.null(keep)) {
        cc <- names(coefs)
        keep <- unlist(lapply(keep, function(x) {
          cc[grepl(x, cc, perl = TRUE, ignore.case=ignore.case)]
        }))
      }
    }
    if (any(unlist(lapply(objects, function(x) is.null(IC(x)))))) {
      ## No iid decomposition/influence functions
      V <- lapply(objects, vcov)
      V <- Reduce(function(...) blockdiag(..., pad=NA), V)
      return(estimate(coef=coefs, vcov=V, keep=keep, ...))
    }
    if (!missing(id) && is.null(id)) { ## Independence between datasets in x,y,...
        nn <- unlist(lapply(
          objects,
          function(x) nrow(x$IC)
        ))
        cnn <- c(0, cumsum(nn))
        id <- list()
        for (i in seq_along(nn)) {
          id <- c(id, list(seq(nn[i]) + cnn[i]))
        }
    }
    if (missing(id)) {
      if (paired) { ## One-to-one dependence between observations in x,y,...
        id <- lapply(objects, function(x) {
          seq_len(NROW(x$IC))
        })
        } else {
            id <- lapply(objects, function(x) x$id)
        }
    } else {
        nn <- unlist(lapply(objects,function(x) NROW(IC(x))))
        if (length(id)==1 && is.logical(id)) {
            if (id) {
                if (any(nn[1]!=nn)) stop("Expected objects of the same size: ", paste(nn,collapse=","))
                id0 <- seq(nn[1]); id <- c()
                for (i in seq(length(nn))) id <- c(id,list(id0))
            } else {
                id <- c()
                N <- cumsum(c(0,nn))
                for (i in seq(length(nn))) id <- c(id,list(seq(nn[i])+N[i]))
            }
        }
        if (length(id)!=length(objects)) stop("Same number of id-elements as model objects expected")
        idlen <- unlist(lapply(id,length))
        if (!identical(idlen,nn)) stop("Wrong lengths of 'id': ", paste(idlen,collapse=","), "; ", paste(nn,collapse=","))
    }
    ids <- ic_all <- c(); count <- 0
    for (z in objects) {
        count <- count+1
        clidx <- NULL
        id0 <- id[[count]]
        icz <- IC(z)
        if (is.null(id0)) {
            id0 <- rownames(icz)
            if (is.null(id0)) stop("Need id for object number ", count)
        }
        if (!missing(subset)) icz <- icz[,subset,drop=FALSE]
        if (!lava.options()$cluster.index) {
            ic0 <- matrix(unlist(by(icz,id0,colSums)),byrow=TRUE,ncol=ncol(icz))
            ids <- c(ids, list(sort(unique(id0))))
        } else {
            if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
            clidx <- mets::cluster.index(id0,mat=icz,return.all=TRUE)
            ic0 <- clidx$X
            ids <- c(ids, list(id0[as.vector(clidx$firstclustid)+1]))
        }
        ic0 <- ic0*NROW(ic0)/length(id0)
        ic_all <- c(ic_all, list(ic0))
    }
    id <- unique(unlist(ids))
    ic0 <- matrix(NA, nrow=length(id),ncol=length(coefs))
    model.index <- c()
    colpos <- 0
    for (i in seq(length(objects))) {
        relpos <- seq_along(coef(objects[[i]]))        
        if (!missing(subset)) relpos <- seq_along(subset)
        ic0[match(ids[[i]], id), relpos + colpos] <- ic_all[[i]]
        midx <- objects[[i]]$model.index
        if (!is.null(midx)) {
          midx <- lapply(midx, function(x) {
            intersect(x, relpos) + colpos
          })
        } else {
          midx <- list(relpos + colpos)
        }
        model.index <- c(model.index, midx)
        colpos <- colpos+tail(relpos,1)
    }
    rownames(ic0) <- id
    ## Rescale each column according to I(obs)/pr(obs)
    for (i in seq(NCOL(ic0))) {
      pr <- mean(!is.na(ic0[,i]))
      ic0[,i] <- ic0[,i]/pr
    }
    ic0[is.na(ic0)] <- 0

    res <- estimate.default(NULL,
      coef = coefs, stack = FALSE, data = NULL,
      IC = ic0, id = id, keep = keep
      )
    if (is.null(keep))
      res$model.index <- model.index
    return(res)
}


##' @export
merge.lm <- function(x,y,...) {
    args <- c(list(x,y),list(...))
    nn <- names(formals(merge.estimate)[-seq(3)])
    idx <- na.omit(match(nn,names(args)))
    models <- args; models[idx] <- NULL
    mm <- lapply(args,function(x) tryCatch(estimate(x),error=function(e) NULL))
    names(mm)[1:2] <- c("x","y")
    ii <- which(unlist(lapply(mm,is.null)))
    if (length(ii)>0) mm[ii] <- NULL
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

