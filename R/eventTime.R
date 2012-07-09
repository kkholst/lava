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
##' regression(m) <- X1
##' regression(m) <- X2
##' @export
##' @param object 
##' @param formula 
##' @param eventName 
##' @param ... 
eventTime <- function(object,formula,eventName,...){
  ff <- as.character(formula)
  timeName <- all.vars(update.formula(formula,"~1"))
  if (length(timeName)==0){
    timeName <- "observedTime"
    rhs <- ff[[2]]
  }else{
    rhs <- ff[[3]]
  }
  rhs <- tolower(rhs)
  latentTimes <- strsplit(rhs,"[(,)]")[[1]]
  if (latentTimes[1]!="min")
    stop(paste("Formula ",formula," does not have the required form, e.g. ~min(T1=1,T2=2,C=0), see (examples in) help(eventTime)."))
  browser()
  latentTimes <- latentTimes[-1]
  NT <- length(latentTimes)
  events <- vector(NT,mode="character")
  for (lt in 1:NT){
    tmp <- strsplit(latentTimes[lt],"=")[[1]]
    stopifnot(length(tmp) %in% c(1,2))
    if (length(tmp)==1){
      events[lt] <- paste("cause",lt,sep=".")
      latentTimes[lt] <- tmp
    }
    else{
      events[lt] <- tmp[2]
      latentTimes[lt] <- tmp[1]
    }
  }
  m <- regression(m,paste("~",timeName))
  m <- regression(m,paste("~",eventName))
  distribution(m,paste("~",timeName)) <- time.lvm(times=latentTimes)
  distribution(m,paste("~",timeName)) <- event.lvm(events=events,times=latentTimes)
  m
}


##' @export
time.lvm <- function(times,...) {
  f <- function(n,times,...) pmin(times)
  return(f)
}
##' @export
event.lvm <- function(events,times,...) {
  f <- function(n,vars) {
    e <- rep(events[1],n)
    pmin(events,times)
  }
  return(f)
}
