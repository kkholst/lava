##' @export
categorical <- function(x,formula,K,beta,p,...) {
    if (missing(beta)) beta <- rep(0,K)
    fname <- paste(gsub(" ","",deparse(formula)),seq(K)-1,sep=":")
    fpar <- names(beta)
    if (is.null(fpar)) fpar <- fname
    y <- getoutcome(formula)
    if (missing(p)) p <- rep(1/K,K-1)
    pname <- names(p)
    if (is.null(pname)) pname <- rep(NA,K-1)
    ordinal(x,K=K,liability=FALSE,p=p,constrain=pname) <- attr(y,"x")
    parameter(x,fpar,start=beta) <- fname
    val <- paste("function(x,p,...) p[\"",fpar[1],"\"]*(x==0)",sep="")
    for (i in seq(K-1)) {        
        val <- paste(val,"+p[\"",fpar[i+1],"\"]*(x==",i,")",sep="")
    }    
    functional(x,formula) <- eval(parse(text=val))
    return(x)
}

##' @export
'categorical<-' <- function(x,...,value) categorical(x,value,...)
