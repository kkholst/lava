##' @export
`%++%.lvm` <- function(x,y) merge(x,y)

##' @export
"+.lvm" <- function(x,...) {
  merge(x,...)
}

## ##' @export
## "+.lm" <- function(x,...) {
##   merge(x,...)
## }

##' @export
merge.lvm <- function(x,y,...) {
  objects <- list(x,y,...)
  if (length(objects)<2) return(x)
  m <- objects[[1]]
  for (i in seq(2,length(objects))) {
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
"+.estimate" <- function(x,...) {
  merge(x, ..., paired=TRUE)
}

##' @export
merge.estimate <- function(x,y,...,id,paired=FALSE,labels=NULL,keep=NULL,subset=NULL) {
    objects <- list(x, estimate(y), ...)
    if (length(nai <- names(objects)=="NA")>0)
    names(objects)[which(nai)] <- ""
    if (!missing(subset)) {
        coefs <- unlist(lapply(objects, function(x) coef(x)[subset]))
    } else {
        coefs <- unlist(lapply(objects,coef))
    }
    if (!is.null(labels)) {
        names(coefs) <- labels
    } else {
        names(coefs) <- make.unique(names(coefs))
    }
    if (!missing(id) && is.null(id)) { ## Independence between datasets in x,y,...
        nn <- unlist(lapply(objects,function(x) nrow(x$iid)))
        cnn <- c(0,cumsum(nn))
        id <- list()
        for (i in seq_along(nn)) id <- c(id,list(seq(nn[i])+cnn[i]))
    }
    if (missing(id)) {
      if (paired) { ## One-to-one dependence between observations in x,y,...
            id <- rep(list(seq(nrow(x$iid))),length(objects))
        } else {
            id <- lapply(objects,function(x) x$id)
        }
    } else {
        nn <- unlist(lapply(objects,function(x) NROW(iid(x))))
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
    ##if (any(unlist(lapply(id,is.null)))) stop("Id needed for each model object")
    ##iid <- Reduce("cbind",lapply(objects,iid))
    ids <- iidall <- c(); count <- 0
    for (z in objects) {
        count <- count+1
        clidx <- NULL
        id0 <- id[[count]]
        iidz <- iid(z)
        if (is.null(id0)) {
            id0 <- rownames(iidz)
            if (is.null(id0)) stop("Need id for object number ", count)
        }
        if (!missing(subset)) iidz <- iidz[,subset,drop=FALSE]
        if (!lava.options()$cluster.index) {
            iid0 <- matrix(unlist(by(iidz,id0,colSums)),byrow=TRUE,ncol=ncol(iidz))
            ids <- c(ids, list(sort(unique(id0))))

        } else {
            if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")
            clidx <- mets::cluster.index(id0,mat=iidz,return.all=TRUE)
            iid0 <- clidx$X
            ids <- c(ids, list(id0[as.vector(clidx$firstclustid)+1]))
        }
        iidall <- c(iidall, list(iid0))
    }
    id <- unique(unlist(ids))
    iid0 <- matrix(NA, nrow=length(id),ncol=length(coefs))
    model.index <- c()
    colpos <- 0
    for (i in seq(length(objects))) {
        relpos <- seq_along(coef(objects[[i]]))        
        if (!missing(subset)) relpos <- seq_along(subset)
        iid0[match(ids[[i]],id),relpos+colpos] <- iidall[[i]]
        model.index <- c(model.index,list(relpos+colpos))
        colpos <- colpos+tail(relpos,1)
    }
    rownames(iid0) <- id
    ## Rescale each column according to I(obs)/pr(obs)
    for (i in seq(NCOL(iid0))) {
      pr <- mean(!is.na(iid0[,i]))
      iid0[,i] <- iid0[,i]/pr
    }
    iid0[is.na(iid0)] <- 0

    res <- estimate.default(NULL, coef=coefs, stack=FALSE, data=NULL,
                            iid=iid0, id=id, keep=keep)
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

