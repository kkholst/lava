##' @export
"ordinal<-" <- function(x,...,value) UseMethod("ordinal<-")

##' @S3method ordinal<- lvm
"ordinal<-.lvm" <- function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(ordinal(x,all.vars(value),...))
  }
  ordinal(x, value, ...)
}

##' @export
"ordinal" <- function(x,...) UseMethod("ordinal")

##' @S3method print ordinal.lvm
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

##' @S3method ordinal lvm
`ordinal.lvm` <- function(x,var=NULL,K=2, constrain, start, ...) {
    if (is.null(var)) {
        ordidx <- unlist(x$attributes$ordinal)
        KK <- unlist(x$attributes$nordinal)
        idx <- x$attributes$ordinalparname
        fix <- lapply(idx,function(z) x$exfix[z])
        if (length(ordidx)>0) {
            val <- names(ordidx)
            return(structure(val,K=KK,idx=idx,fix=fix,class="ordinal.lvm"))
        }
        else
            return(NULL)
    }
    if (length(var)>length(K)) K <- rep(K[1],length(var))
    if (length(var)==1 && !missing(constrain)) constrain <- list(constrain)
    for (i in seq_len(length(var))) {
        if (K[i]>2) {
            parname <- paste(var[i],":",paste(seq(K[i]-1)-1,seq(K[i]-1),sep="|"),sep="")
            newpar <- if (missing(start))
                rep(-1,K[i]-1) else start[[i]]
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
        x$attributes$type[var[i]] <- ifelse(K[i]>2,"Ordinal","Binary")
        if (K[i]>2) intfix(x,var[i],NULL) <- 0
    }
    x$attributes$ordinal[var] <- TRUE
    x$attributes$nordinal[var] <- K
    x$attributes$normal[var] <- FALSE
    covfix(x,var,NULL) <- 1
  return(x)
}
