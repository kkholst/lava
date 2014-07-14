color.ordinal <- function(x,subset=vars(x),...) {
    return(list(vars=intersect(subset,ordinal(x)),col="indianred1"))
}

ordinal.sim.hook <- function(x,data,p,...) {
    ovar <- ordinal(x)
    for (i in seq_len(length(ovar))) {
        if (attributes(ovar)$liability[i]) {
            idx <- attributes(ovar)$idx[[ovar[i]]]
            if (length(idx)==0) {
                breaks <- c(-Inf,0,Inf)
            } else {
                breaks <- c(-Inf,ordreg_threshold(p[idx]),Inf)
            }
            z <- cut(data[,ovar[i]],breaks=breaks)
            data[,ovar[i]] <- as.numeric(z)-1
        }
        K <- attributes(ovar)$K[i]
        lab <- attributes(ovar)$labels[ovar[i]][[1]]
        if (!is.null(lab))
            data[,ovar[i]] <- factor(data[,ovar[i]],
                                     levels=seq(K)-1,
                                     labels=lab)
        
    }
    return(data)
}
##addhook("ordinal.sim.hook","sim.hooks")

##' @export
"ordinal<-" <- function(x,...,value) UseMethod("ordinal<-")

##' @export
"ordinal<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(ordinal(x,all.vars(value),...))
  }
  ordinal(x, value, ...)
}

##' @export
"ordinal" <- function(x,...) UseMethod("ordinal")

##' @export
print.ordinal.lvm <- function(x,...) {
  cat(rep("_",28),"\n",sep="")
  for (i in x) {
    val <- attr(x,"fix")[[i]]
    if (length(val)==0)
      cat(paste(i,"binary",sep=":"),"\n")
    else print(unlist(attr(x,"fix")[[i]]),quote=FALSE)
    cat(rep("_",28),"\n",sep="")
  }
}

##' @export
`ordinal.lvm` <- function(x,var=NULL,K=2, constrain, breaks=NULL, p, liability=TRUE, labels, ...) {
    if (is.null(var)) {
        ordidx <- unlist(x$attributes$ordinal)
        KK <- unlist(x$attributes$nordinal)
        idx <- x$attributes$ordinalparname
        fix <- lapply(idx,function(z) x$exfix[z])
        liability <- x$attributes$liability
        labels <- x$attributes$labels
        if (length(ordidx)>0) {
            val <- names(ordidx)
            return(structure(val,K=KK,idx=idx,fix=fix,liability=liability,labels=labels,class="ordinal.lvm"))
        }
        else
            return(NULL)
    }
    if (!missing(p)) breaks <- qnorm(cumsum(p))
    if (!is.null(breaks)) {
        breaks <- ordreg_ithreshold(breaks)
        K <- length(breaks)+1
    }
    if (length(var)>length(K)) K <- rep(K[1],length(var))
    if (length(var)==1 && !missing(constrain)) constrain <- list(constrain)

    
    addvar(x) <- var
    for (i in seq_len(length(var))) {
        if (K[i]>2) {
            parname <- paste(var[i],":",paste(seq(K[i]-1)-1,seq(K[i]-1),sep="|"),sep="")
            newpar <- if (is.null(breaks)) {
                rep(-1,K[i]-1)
            } else if (is.list(breaks)) breaks[[i]] else breaks
            if (length(newpar)<K[i]-1) stop("Wrong number of starting values")
            newfix <- if (missing(constrain))
                rep(list(NA),length(newpar)) else constrain[[i]]
            if (length(newfix)<K[i]-1) stop("Wrong number of constraints")
            names(newpar) <- names(newfix) <- parname
            
            parameter(x,newfix,start=newpar) <- names(newfix)
            ## pex <- parname%in%names(x$expar)
            ## if (any(pex <- parname%in%names(x$expar))) {
            ##     if (!all(pex)) stop("Cannot change number of categories! Re-specify model.")
            ##     x$attributes$iordinal[var] <- list(idx)       
            ## }
            x$attributes$ordinalparname[var[i]] <- list(names(newfix))
        }
        x$attributes$type[var[i]] <- ifelse(K[i]>2,"ordinal","binary")
        if (K[i]>2) intfix(x,var[i],NULL) <- 0
        if (!liability) {
            mytr <- function(y,p,idx,...) {
                breaks <- c(-Inf,ordreg_threshold(p[idx]),Inf)
                as.numeric(cut(y,breaks=breaks))-1
            }
            myalist <- substitute(alist(y=,p=,idx=pp),
                                  list(pp=x$attributes$ordinalparname[[var[i]]]))
            formals(mytr) <- eval(myalist)
            transform(x,var[i],post=FALSE) <- mytr
            
        }
    }
    x$attributes$liability[var] <- liability
    x$attributes$ordinal[var] <- TRUE
    if (!missing(labels)) x$attributes$labels[var] <- list(labels)
    x$attributes$nordinal[var] <- K
    x$attributes$normal[var] <- FALSE
    covfix(x,var,NULL) <- 1
    if (is.null(index(x))) index(x) <- reindex(x)
    return(x)
}
