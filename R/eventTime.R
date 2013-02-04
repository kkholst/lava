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
eventTime <- function(object,formula,eventName,...){
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
  for (lt in 1:NT){
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
  addvar(object) <- timeName
  ## m <- regression(m,formula(paste("~",timeName,sep="")))
  if (missing(eventName)) eventName <- "Event"
  eventTime <- list(names=c(timeName,eventName),
                    latentTimes=gsub(" ","",latentTimes),
                    events=gsub(" ","",events))
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
coxWeibull.lvm <- function(shape=1,scale) {
  f <- function(n,mu,var,...) {
    (- (log(runif(n)) * (1 / scale) * exp(-mu)))^(1/shape)
  }
  attr(f,"family") <- list(family="weibull",par=c(shape,scale))
  return(f)
}

##' @export
coxExponential.lvm <- function(scale){
  coxWeibull.lvm(shape=1,scale)
}

##' @export
coxGompertz.lvm <- function(shape=1,scale) {
  f <- function(n,mu,var,...) {
    (1/shape) * log(1 - (shape/scale) * (log(runif(n)) * exp(-mu)))
  }
  attr(f,"family") <- list(family="gompertz",par=c(shape,scale))
  return(f)
}



