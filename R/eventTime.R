##' Add an observed event time outcome to a latent variable model.
##'
##' For example, if the model 'm' includes latent event time variables
##' are called 'T1' and 'T2' and 'C' is the end of follow-up (right censored), 
##' then one can specify 
##'
##' \code{eventTime(object=m,formula=ObsTime~min(T1=a,T2=b,C=0,"ObsEvent"))}
##'
##' when data are simulated from the model
##' one gets 2 new columns: 
##'
##' - "ObsTime": the smallest of T1, T2 and C
##' - "ObsEvent": 'a' if T1 is smallest, 'b' if T2 is smallest and '0' if C is smallest
##'
##' Note that "ObsEvent" and "ObsTime" are names specified by the user.
##' 
##' @author Thomas A. Gerds
##' @keywords survival models regression
##' @examples
##'
##' m <- lvm() 
##' distribution(m,~X2) <- binomial.lvm()
##' regression(m) <- T1~f(X1,-.5)+f(X2,0.3)
##' regression(m) <- T2~f(X2,0.6)
##' distribution(m,~T1) <- coxWeibull.lvm(scale=1/100)
##' distribution(m,~T2) <- coxWeibull.lvm(scale=1/100)
##' distribution(m,~C) <- coxWeibull.lvm(scale=1/100)
##' m <- eventTime(m,time~min(T1=1,T2=2,C=0),"event")
##' d <- sim(m,10000)
##'
##' 
##' m <- eventTime(m,cens.otime~min(T1,T2=E2,C=0),"event")
##' sim(m,10)
##' @export
##' @param object Model object
##' @param formula Formula (see details)
##' @param eventName Event names
##' @param \dots Additional arguments to lower levels functions
eventTime <- function(object,formula,eventName,...) {
  if (missing(formula)) return(object$attributes$eventHistory)
  ff <- as.character(formula)
  timeName <- all.vars(update.formula(formula,"~1"))
  if (length(timeName)==0){
    timeName <- "observedTime"
    rhs <- ff[[2]]
  }else{
    rhs <- ff[[3]]
  }
  ## rhs <- tolower(rhs)
  latentTimes <- strsplit(rhs,"[(,)]")[[1]]
  if (latentTimes[1]!="min")
    stop(paste("Formula ",formula," does not have the required form, e.g. ~min(T1=1,T2=2,C=0), see (examples in) help(eventTime)."))
  latentTimes <- latentTimes[-1]
  NT <- length(latentTimes)
  events <- vector(NT,mode="character")
  for (lt in seq_len(NT)){
    tmp <- strsplit(latentTimes[lt],"=")[[1]]
    stopifnot(length(tmp) %in% c(1,2))
    if (length(tmp)==1){
      events[lt] <- as.character(lt)
      latentTimes[lt] <- tmp
    }
    else{
      events[lt] <- tmp[2]
      latentTimes[lt] <- tmp[1]
    }
  }
  events <- gsub(" ","",events)
  suppressWarnings(eventnum <- as.numeric(events)) 
  if (all(!is.na(eventnum))) {
    events <- eventnum
  } else {
    events <- gsub("\"","",events)
  }

  addvar(object) <- timeName
  ## m <- regression(m,formula(paste("~",timeName,sep="")))
  if (missing(eventName)) eventName <- "Event"
  eventTime <- list(names=c(timeName,eventName),
                    latentTimes=gsub(" ","",latentTimes),
                    events=events
                    )
  if (is.null(object$attributes$eventHistory)) {
    object$attributes$eventHistory <- list(eventTime)
    names(object$attributes$eventHistory) <- timeName
  } else {
    object$attributes$eventHistory[[timeName]] <- eventTime
  }
  return(object)
}


## addhook("color.eventHistory","color.hooks")
## color.eventHistory <- function(x,subset=vars(x),...) {
##   return(list(vars=intersect(subset,binary(x)),col="indianred1"))
## }

addhook("plothook.eventHistory","plot.post.hooks")
plothook.eventHistory <- function(x,...) {
  eh <- x$attributes$eventHistory
  ehnames <- unlist(lapply(eh,function(x) x$names))
  for (f in eh) {
    x <- regression(x,to=f$names[1],from=f$latentTimes)
    latent(x) <- f$latentTimes
  }
  return(x)
}


addhook("print.eventHistory","print.hooks")
print.eventHistory <- function(x,...) { 
  if (is.null(eh <- x$attributes$eventHistory)) return(NULL)
  ehnames <- unlist(lapply(eh,function(x) x$names))
  cat("Event History Model\n")
  ff <- formula(x,TRUE)
  R <- c()
  for (f in ff) {
    oneline <- as.character(f);
    y <- gsub(" ","",strsplit(f,"~")[[1]][1])
    if (!(y %in% ehnames)) {
      col1 <- as.character(oneline)
      D <- attributes(distribution(x)[[y]])$family
      col2 <- "Normal"
      if (!is.null(D$family)) col2 <- paste(D$family,sep="")
      if (!is.null(D$link)) col2 <- paste(col2,"(",D$link,")",sep="")
      if (!is.null(D$par)) col2 <- paste(col2,"(",paste(D$par,collapse=","),")",sep="")      
      R <- rbind(R,c(col1,"  ",col2))
    }
  }
  for (y in names(eh)) {
    col1 <- paste(y, " = min(",paste(eh[[y]]$latentTimes,collapse=","),")",sep="")
    eh[[y]]$names[2]
    col2 <- paste(eh[[y]]$names[2], " := {",paste(eh[[y]]$events,collapse=","),"}",sep="")
    R <- rbind(R,c(col1,"",col2))
  }
  rownames(R) <- rep("",nrow(R)); colnames(R) <- rep("",ncol(R))
  print(R,quote=FALSE,...)
  cat("\n")
  TRUE
}

addhook("simulate.eventHistory","sim.hooks")

simulate.eventHistory <- function(x,data,...){
  if (is.null(eventTime(x))) {
    return(data)
  }
  else{
    for (eh in eventTime(x)) {
      if (any((found <- match(eh$latentTimes,names(data),nomatch=0))==0)){
        warning("Cannot find latent time variable: ",
                eh$latentTimes[found==0],".")
      }
      else{
        for (v in 1:length(eh$latentTimes)) {
          if (v==1){ ## initialize with the first latent time and event
            eh.time <- data[,eh$latentTimes[v]]
            eh.event <- rep(eh$events[v],NROW(data))
          } else{ ## now replace if next time is smaller
            ## in case of tie keep the first event
            eh.event[data[,eh$latentTimes[v]]<eh.time] <- eh$events[v]
            eh.time <- pmin(eh.time,data[,eh$latentTimes[v]])
          }
        }
      }
      data[,eh$names[1]] <- eh.time
      data[,eh$names[2]] <- eh.event
    }
    return(data)
  }
}



##' @export
coxWeibull.lvm <- function(shape=1,scale,rate=1/scale) {
  f <- function(n,mu,var,...) {
    (- (log(runif(n)) * rate * exp(-mu)))^(1/shape)
  }
  attr(f,"family") <- list(family="weibull",par=c(shape,scale))
  return(f)
}


##' Add time-varying covariate effects to model
##' 
##' @title Time-dependent parameters
##' @param object Model
##' @param formula Formula with rhs specifying time-varying covariates
##' @param rate Optional (log)-rate(ratio) parameters
##' @param timecut Time intervals
##' @param type Type of model (default piecewise constant intensity)
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
##' @examples
##' 
##' ## Piecewise constant hazard
##' m <- lvm(y~1)
##' m <- timedep(m,y~1,timecut=c(0,5),rate=c(0.5,0.3))
##' 
##' \dontrun{
##' d <- sim(m,1e4); d$status <- TRUE
##' dd <- mets::lifetable(Surv(y,status)~1,data=d,breaks=c(5));
##' exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval, dd, family=poisson)))
##' }
##' 
##' 
##' ## Piecewise constant hazard and time-varying effect of z1
##' m <- lvm(y~1)
##' distribution(m,~z1) <- ones.lvm(0.5)
##' R <- log(cbind(c(0.2,0.7,0.9),c(0.5,0.3,0.3)))
##' m <- timedep(m,y~z1,timecut=c(0,3,5),rate=R)
##' 
##' \dontrun{
##' d <- sim(m,1e4); d$status <- TRUE
##' dd <- lifetable(Surv(y,status)~z1,data=d,breaks=c(3,5));
##' exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval+z1:interval, dd, family=poisson)))
##' }
##' 
##' 
##' 
##' ## Explicit simulation of time-varying effects
##' m <- lvm(y~1)
##' distribution(m,~z1) <- ones.lvm(0.5)
##' distribution(m,~z2) <- binomial.lvm(p=0.5)
##' #variance(m,~m1+m2) <- 0
##' #regression(m,m1[m1:0] ~ z1) <- log(0.5)
##' #regression(m,m2[m2:0] ~ z1) <- log(0.3)
##' regression(m,m1 ~ z1,variance=0) <- log(0.5)
##' regression(m,m2 ~ z1,variance=0) <- log(0.3)
##' intercept(m,~m1+m2) <- c(-0.5,0)
##' m <- timedep(m,y~m1+m2,timecut=c(0,5))
##' 
##' \dontrun{
##' d <- sim(m,1e5); d$status <- TRUE
##' dd <- lifetable(Surv(y,status)~z1,data=d,breaks=c(5))
##' exp(coef(glm(events ~ offset(log(atrisk)) + -1 + interval + interval:z1, dd, family=poisson)))
##' }
timedep <- function(object,formula,rate,timecut,type="coxExponential.lvm",...) {
    if (missing(timecut)) stop("'timecut' needed")
    ##if (inherits(formula,"formula"))
    ff <- getoutcome(formula)
    simvars <- attributes(ff)$x
    if (is.null(object$attributes$simvar)) {
        object$attributes$simvar <- list(simvars)
        names(object$attributes$simvar) <- ff
    } else {
        object$attributes$simvar[[ff]] <- simvars
    }
    if (missing(rate)) rate <- rep(1,length(timecut))
    args <- list(timecut=timecut,rate=rate,...)
    distribution(object,ff) <- do.call(type,args)
    return(object)
}


##' @export
coxExponential.lvm <- function(scale=1,rate,timecut){
    if (missing(rate)) rate=1/scale
    if (missing(scale)) scale=1/rate
    if (missing(timecut)) {
        return(coxWeibull.lvm(shape=1,scale))
    }

    if (NROW(rate)>length(timecut))
        stop("Number of time-intervals (cuts) does not agree with number of rate parameters (beta0)")    
    par <- paste(timecut,rate,sep=":")
    if (is.matrix(rate)) par <- "..."
    timecut <- c(timecut,Inf)
    f <- function(n,mu,...) {
        Ai <- function() {
            vals <- matrix(0,ncol=length(timecut)-1,nrow=n)
            ival <- numeric(n)
            if (is.matrix(rate)) {
                mu <- cbind(mu[,1],cbind(1,as.matrix(mu[,-1]))%*%t(rate))
                rate <- rep(1,length(timecut)-1)                
            }   
            for (i in seq(length(timecut)-1)) {
                u <- -log(runif(n)) ##rexp(n,1)
                if (NCOL(mu)>1) {
                    vals[,i] <-  timecut[i] + u*exp(-mu[,1]-mu[,i+1])/(rate[i])
                } else {
                    vals[,i] <-  timecut[i] + u*exp(-mu)/(rate[i])
                }
                idx <- which(vals[,i]<=timecut[i+1] & ival==0)
                ival[idx] <- vals[idx,i]
            }
            ival
        }
        Ai()
    }
    attributes(f)$family <- list(family="CoxExponential",par=par)
    return(f)
}

##' @export
aalenExponential.lvm <- function(scale=1,rate,timecut=0){
    if (missing(rate)) rate=1/scale
    if (missing(scale)) scale=1/rate
    if (missing(timecut)==1) {
        return(coxWeibull.lvm(shape=1,scale))
    }
    
    if (length(rate)>length(timecut))
        stop("Number of time-intervals (cuts) does not agree with number of rate parameters (beta0)")
    par <- paste(timecut,rate,sep=":")
    if (is.matrix(rate)) par <- "..."
    timecut <- c(timecut,Inf)
    f <- function(n,mu,...) {
        Ai <- function() {
            vals <- matrix(0,ncol=length(timecut)-1,nrow=n)
            ival <- numeric(n)
            if (is.matrix(rate)) {
                mu <- cbind(mu[,1],cbind(1,as.matrix(mu[,-1]))%*%t(rate))
                rate <- rep(1,length(timecut)-1)                
            }
            for (i in seq(length(timecut)-1)) {
                u <- -log(runif(n)) ##rexp(n,1)
                if (NCOL(mu)>1) {
                    vals[,i] <-  timecut[i] + u/(rate[i]+mu[,1]+mu[,i+1])
                } else {
                    vals[,i] <-  timecut[i] + u/(rate[i]+mu)
                }
                idx <- which(vals[,i]<=timecut[i+1] & ival==0)
                ival[idx] <- vals[idx,i]
            }
            ival
        }
        Ai()
    }
    attributes(f)$family <- list(family="aalenExponential",par=par)
    return(f)
}


##' @export
coxGompertz.lvm <- function(shape=1,scale) {
  f <- function(n,mu,var,...) {
    (1/shape) * log(1 - (shape/scale) * (log(runif(n)) * exp(-mu)))
  }
  attr(f,"family") <- list(family="gompertz",par=c(shape,scale))
  return(f)
}



