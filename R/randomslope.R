##' @export
"randomslope<-" <- function(x,...,value) UseMethod("randomslope<-")

##' @S3method randomslope<- lvm
"randomslope<-.lvm" <- function(x, ..., value) {
  randomslope(x, covar=value, ...)
}

##' @export
`randomslope` <-
function(x,...) UseMethod("randomslope")

##' @S3method randomslope lvm
`randomslope.lvm` <-
function(x,covar,random=NULL,response=NULL,param,postfix,clear=FALSE,zero=TRUE,...) {
  if (missing(covar)) {
    rsidx <- unlist(nodeData(Graph(x), attr="randomslope"))
    if (length(rsidx)>0)
      return(names(rsidx)[rsidx])
    else
      NULL    
  }
  if (class(covar)[1]=="formula") {
    covar <- all.vars(covar)
  }
  if (clear) {
    ##    x <- addattr(x,attr="shape",var=var,val="rectangle")
    nodeData(Graph(x), covar, attr="randomslope") <- FALSE
  } else {
    if (!is.null(random) & !is.null(response)) {
      if (class(random)[1]=="formula") {
        random <- all.vars(random)
      }
      if (class(response)[1]=="formula") {
        response <- all.vars(response)
      }
      if (length(covar)!=length(response)) stop("Vectors should be of the same length!")
      if (!(random%in%latent(x))) {
        addvar(x) <- random
        latent(x) <- random
      }
      if (missing(param) || !is.null(param)) {
        if (!missing(postfix))
          newlatent <-  paste(random,postfix,sep="")
        else
          newlatent <-  paste(random,covar,sep=".")
        covariance(x,random) <- 1
        for (i in 1:length(covar)) {
          if (missing(param)) {
            x <- regression(x,to=newlatent[i],from=random)
          } else {
            if (class(param)[1]=="formula") {
              param <- all.vars(param)
            }    
            if (length(param)!=length(newlatent))
              param <- rep(param,length(newlatent))
            regfix(x,to=newlatent[i], from=random) <- param[i]
          }              
          regfix(x,to=response[i],from=newlatent[i]) <- covar[i]
          latent(x) <- newlatent[i]
          covariance(x,newlatent[i]) <- 0
        }
      } else {
        for (i in 1:length(covar)) {
          regfix(x,to=response[i],from=random) <- covar[i]
        }
      }
    } else {
      nodeData(Graph(x), covar, attr="randomslope") <- TRUE
    }
  }
  index(x) <- reindex(x)
  return(x)
}

##' @S3method randomslope lvmfit
`randomslope.lvmfit` <-
  function(x,...) {
    randomslope(Model(x),...)
  }

