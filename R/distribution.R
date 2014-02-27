###{{{ distribution

##' @export
"distribution<-" <- function(x,...,value) UseMethod("distribution<-")

##' @export
"distribution" <- function(x,...,value) UseMethod("distribution")

##' @S3method distribution<- lvm
"distribution<-.lvm" <- function(x,variable,...,value) {
  if (class(variable)[1]=="formula")
    variable <- all.vars(variable)  
  if (length(variable)==1) {
    addvar(x) <- as.formula(paste("~",variable))
    if (is.numeric(value)) value <- list(value)
    x$attributes$distribution[[variable]] <- value ##eval(parse(text=mytext))
    return(x)
  }    
  if ((length(value)!=length(variable) & length(value)!=1))
    stop("Wrong number of values")
  for (i in 1:length(variable))
    if (length(value)==1) {
      distribution(x,variable[i],...) <- value
    } else {
      distribution(x,variable[i],...) <- value[[i]]
    }
  return(x)
  
}

##' @S3method distribution lvm
"distribution.lvm" <- function(x,var,...) {
  x$attributes$distribution[var]
}

###}}} distribution

###{{{ normal/gaussian

##' @export
normal.lvm <- function(link="identity",mean,sd,log=FALSE,...) {
  rnormal <- if(log) rlnorm else rnorm
  fam <- gaussian(link); fam$link <- link
  if (!missing(mean) & !missing(sd)) 
    f <- function(n,mu,var,...) rnormal(n,fam$linkinv(mean),sd)
  else
    f <- function(n,mu,var,...) {      
      rnormal(n,fam$linkinv(mu),sqrt(var))
    }
  attr(f,"family") <- fam
  return(f)  
}

##' @export
gaussian.lvm <- normal.lvm

##' @export
lognormal.lvm <- function(...) normal.lvm(...,log=TRUE)


###}}} normal/gaussian

###{{{ poisson

##' @export
poisson.lvm <- function(link="log",lambda,...) {  
  fam <- poisson(link); fam$link <- link
 if (!missing(lambda))
    f <- function(n,mu,...) rpois(n,lambda)
 else
   f <- function(n,mu,...) {
     if (missing(n)) {
       return(fam)
     }
     rpois(n,fam$linkinv(mu))
   }
  attr(f,"family") <- fam
  attr(f,"var") <- FALSE
  return(f)  
} 

###}}} poisson

###{{{ threshold

threshold.lvm <- function(p,labels=NULL,...) {
    if (sum(p)>1 || any(p<0 | p>1)) stop("wrong probability vector") ;
    if (!missing(labels))
    return(function(n,...) {
        return(cut(rnorm(n),breaks=c(-Inf,qnorm(cumsum(p)),Inf),labels=labels))
    })
    function(n,...)
        cut(rnorm(n),breaks=c(-Inf,qnorm(cumsum(p)),Inf))    
}

###}}}

###{{{ binomial

##' @export
binomial.lvm <- function(link="logit",p,size=1) {
  fam <- binomial(link); fam$link <- link
  if (!missing(p))
    f <- function(n,mu,var,...) rbinom(n,size,p)
  else {
    f <- function(n,mu,var,...) {
      if (missing(n)) {
        return(fam)
      }
      rbinom(n,size,fam$linkinv(mu))
    }
    ## f <- switch(link,
    ##             logit = 
    ##             function(n,mu,var,...) rbinom(n,1,tigol(mu)),
    ##             cloglog =
    ##             function(n,mu,var,...) rbinom(n,1,1-exp(-exp(1-mu))),
    ##             function(n,mu,var=1,...) rbinom(n,1,pnorm(mu,sd=sqrt(var)))
    ##             ### function(n,mu=0,var=1,...) (rnorm(n,mu,sqrt(var))>0)*1
    ##             )    
  }
  attr(f,"family") <- fam
  attr(f,"var") <- FALSE
  return(f)
}

##' @export
logit.lvm <- binomial.lvm("logit")

##' @export
probit.lvm <- binomial.lvm("probit")

###}}} poisson

###{{{ Gamma


##' @export
Gamma.lvm <- function(link="inverse",shape,rate,unit=FALSE,var=FALSE,log=FALSE,...) {
  fam <- Gamma(link); fam$link <- link
  rgam <- if (!log) rgamma else function(...) log(rgamma(...))
  if (!missing(shape) & !missing(rate))
    f <- function(n,mu,var,...) rgam(n,shape=shape,rate=rate)
  if (!missing(shape) & missing(rate)) {
    if (unit) 
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=shape)
    else if (var) 
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=sqrt(shape/var))
    else 
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=shape/fam$linkinv(mu))
  }
  if (missing(shape) & !missing(rate)) {
    if (unit) 
      f <- function(n,mu,var,...) rgam(n,shape=shape,rate=rate)
    else if (var) 
      f <- function(n,mu,var,...) rgam(n,shape=rate^2*var,rate=rate)
    else
      f <- function(n,mu,var,...) rgam(n,shape=rate*fam$linkinv(mu),rate=rate)
  }
  if (missing(shape) & missing(rate)) {
    if (var)
      f <- function(n,mu,var,...) rgam(n,shape=var,rate=1)
    else
      f <- function(n,mu,var,...) rgam(n,shape=fam$linkinv(mu),rate=1)
  }
  attr(f,"family") <- fam
  attr(f,"var") <- FALSE
  return(f)  
} 

##' @export
loggamma.lvm <- function(...) Gamma.lvm(...,log=TRUE)


###}}} Gamma

###{{{ student (t-distribution)

##' @export
student.lvm <- function(df=2,mu,sigma,...) {
  if (!missing(mu) & !missing(sigma)) 
    f <- function(n,mu,var,...) mu+sigma*rt(n,df=df)
  else
    f <- function(n,mu,var,...) mu + sqrt(var)*rt(n,df=df)  
  return(f)
}

###}}} student (t-distribution)

###{{{ uniform

##' @export
uniform.lvm <- function(a,b) {
  if (!missing(a) & !missing(b)) 
    f <- function(n,mu,var,...) runif(n,a,b)
  else
    f <- function(n,mu,var,...)
      (mu+(runif(n,-1,1)*sqrt(12)/2*sqrt(var)))
  return(f)
}

###}}} uniform

###{{{ weibull

##' @export
weibull.lvm <- function(scale=1.25,shape=2,cens=Inf,breakties=0) {
  require(survival)
  lambda <- 1/scale
  f <- function(n,mu,var,...) {
    a0 <- function(t) lambda*shape*(lambda*t)^(shape-1)
    A0 <- function(t) (lambda*t)^shape
    A0i <- function(eta) eta^(1/shape)/lambda
    U <- rexp(n, 1) #give everyone a random death time, on the CH scale
    Z <- U*exp(-mu)
    T <- A0i(Z)
    if (breakties!=0)
      T <- T+runif(n,0,breakties)
    if (is.function(cens))
      cens <- cens(n,...)
    if (is.finite(cens[1])) {
      Delta <- (T<cens)
      if (any(!Delta)) {
        T[!Delta] <- cens[!Delta]
      S <- Surv(T,Delta*1)      
      } else {
        S <- T
      }
    } else {
      S <- T
    }
    return(S)
  }
  return(f)
}

###}}} weibull

###{{{ sequence
##' @export
sequence.lvm <- function(a=0,b=1) {
    f <- function(n,...)
      seq(a,b,length.out=n)
  return(f)
}
###}}} sequence

###{{{ ones
##' @export
ones.lvm <- function(fraction=0) {
    f <- function(n,...) {
        val <- rep(1,n)
        if (fraction>0 && fraction<1) val[seq(n*(1-fraction))] <- 0
        val
        }
  return(f)
}
###}}} ones
