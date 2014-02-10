##' @S3method %+% lvm
`%+%.lvm` <- function(x,y) merge(x,y)

##' @S3method merge lvm
merge.lvm <- function(x,y,...) {  
  objects <- list(x,y,...)  
  if (length(objects)<2) return(x)
  m <- objects[[1]]
  for (i in 2:length(objects)) {
    m2 <- objects[[i]]
    if (length(latent(m2))>0)
      latent(m) <- latent(m2)
    if (length(m2$constrain)>0)
      m$constrain <- c(m$constrain,m2$constrain)
    M <- (index(m2)$A)
    P <- (index(m2)$P)
    nn <- vars(m2)
    for (j in 1:nrow(M)) {
      if (any(idx <- M[j,]!=0)) {
        val <- as.list(rep(NA,sum(idx==TRUE)))
        if (any(idx. <- !is.na(m2$par[j,idx])))
          val[idx.] <- m2$par[j,idx][idx.]
        if (any(idx. <- !is.na(m2$fix[j,idx])))
          val[idx.] <- m2$fix[j,idx][idx.]
        regression(m,nn[idx],nn[j],silent=TRUE) <- val
      }
      P0 <- P[j,]; P0[seq_len(j-1)] <- 0
        if (any(idx <- P[j,]!=0)) {
          val <- as.list(rep(NA,sum(idx==TRUE)))
          if (any(idx. <- !is.na(m2$covpar[j,idx])))
            val[idx.] <- m2$covpar[j,idx][idx.]
          if (any(idx. <- !is.na(m2$covfix[j,idx])))
            val[idx.] <- m2$covfix[j,idx][idx.]
          covariance(m,nn[idx],nn[j],silent=TRUE) <- val
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
  
##' @S3method merge lm
merge.lm <- function(x,...) {
    lava:::merge.estimate(x,...)
}

##' ##' @S3method merge glm
merge.glm <- function(x,...) {
    lava:::merge.estimate(x,...)
}

##' @S3method merge estimate
merge.estimate <- function(x,y,...,id) {
    objects <- list(x,y, ...)
    coefs <- unlist(lapply(objects,coef))
    names(coefs) <- make.unique(names(coefs))

    if (missing(id)) {
        id <- lapply(objects,function(x) x$id)
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
    if (any(unlist(lapply(id,is.null)))) stop("Id needed for each model object")
    ##iid <- Reduce("cbind",lapply(objects,iid))    
    ids <- iidall <- c(); count <- 0
    for (z in objects) {
        count <- count+1
        clidx <- NULL
        id0 <- id[[count]]
        if (inherits(try(find.package("mets"),silent=TRUE),"try-error")) {
            iid0 <- matrix(unlist(by(iid(z),id0,colSums)),byrow=TRUE,ncol=length(coef(z)))
            ids <- c(ids, list(sort(unique(id0))))

        } else {
            clidx <- mets::cluster.index(id0,mat=iid(z),return.all=TRUE)
            iid0 <- clidx$X
            ids <- c(ids, list(id0[as.vector(clidx$firstclustid)+1]))
        }
        iidall <- c(iidall, list(iid0))
    }
    id <- unique(unlist(ids))
    iid0 <- matrix(0,nrow=length(id),ncol=length(coefs))
    colpos <- 0
    for (i in seq(length(objects))) {
        relpos <- seq(length(coef(objects[[i]])))
        iid0[match(ids[[i]],id),relpos+colpos] <- iidall[[i]]
        colpos <- colpos+tail(relpos,1)
    }
    rownames(iid0) <- id
    estimate.default(NULL, coef=coefs, stack=FALSE, data=NULL, iid=iid0, id=id)
}


