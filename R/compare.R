comparepair <- function(x1,x2) {
  ##if (class(x1)!="lvmfit" | class(x2)!="lvmfit") stop("'lvmfit' object expected.")
  l1 <- do.call("logLik",list(x1),envir=parent.frame(2))
  l2 <- do.call("logLik",list(x2),envir=parent.frame(2))
  ##l1 <- logLik(x1);  l2 <- logLik(x2)
  df1 <- attributes(l1)$df;  df2 <- attributes(l2)$df;
  ##Q <- -2*ifelse(df1<df2, l1-l2, l2-l1); names(Q) <- "chisq"
  Q <- abs(2*(l1-l2))
  names(Q) <- "chisq"
  df <- abs(df1-df2); names(df) <- "df"
  p <- 1-pchisq(Q,df=df)
  values <- c(l1,l2); names(values) <- c("log likelihood (model 1)", "log likelihood (model 2)")

  res <- list(statistic = Q, parameter = df,
              p.value= p, method = "Likelihood ratio test",
              estimate = values)
  class(res) <- "htest"
  return(res)
}


`compare` <-
  function(object,...) UseMethod("compare")

compare.default <- function(object,...,par,contrast,null,scoretest,Sigma) {
  if (!missing(par)) {
    contrast <- rep(0,length(coef(object)))
    myidx <- parpos(Model(object),p=par)        
    contrast[myidx] <- 1
    contrast <- diag(contrast)[contrast!=0,]
    if (!missing(null) && length(null)>1) null <- null[attributes(myidx)$ord]
  }
  ### Wald test
  if (!missing(contrast)) {
    B <- contrast    
    p <- pars(object)
    if (is.vector(B)) { B <- rbind(B); colnames(B) <- names(contrast) }
    if (missing(Sigma)) {
      Sigma <- vcov(object)
    }
    if (ncol(B)<length(p)) {
      nn <- colnames(B)
      myidx <- parpos(Model(object),p=nn)
      B0 <- matrix(0,nrow=nrow(B),ncol=length(coef(object)))
      B0[,myidx] <- B[,attributes(myidx)$ord]
      B <- B0
    }
    if (missing(null)) null <- 0
    Q <- t(B%*%p-null)%*%Inverse(B%*%Sigma%*%t(B))%*%(B%*%p-null)
    df <- qr(B)$rank; names(df) <- "df"
    attributes(Q) <- NULL; names(Q) <- "chisq";
    pQ <- ifelse(df==0,NA,1-pchisq(Q,df))
    method = "Wald test";
    ##    hypothesis <-
    res <- list(##data.name=hypothesis,
                statistic = Q, parameter = df,
                p.value=pQ, method = method
                )
    class(res) <- "htest"
    attributes(res)$B <- B
    return(res)        
  }

  ### Score test
  if (!missing(scoretest)) {
    altmodel <- Model(object)
    if (class(scoretest)[1]=="formula") scoretest <- list(scoretest)
    for (i in scoretest) {
      regression(altmodel) <- i
    }
    p0 <- numeric(length(coef(altmodel)))        
    idx <-  match(coef(Model(object)),coef(altmodel))
    p0[idx] <- coef(object)
    Sc2 <- score(altmodel,p=p0,data=model.frame(object),weigth=Weight(altmodel),
                 estimator=object$estimator,...)
    I <- information(altmodel,p=p0,n=object$data$n,
                     data=model.frame(object),weigth=Weight(object),
                     estimator=object$estimator,...
                     )
    iI <- try(solve(I), silent=TRUE)
    Q <- ifelse (inherits(iI, "try-error"), NA, ## Score test
                 ## rbind(Sc)%*%iI%*%cbind(Sc)
                 (Sc2)%*%iI%*%t(Sc2)
                 )
    attributes(Q) <- NULL; names(Q) <- "chisq"
    df <- length(p0)-length(coef(object)); names(df) <- "df"
    pQ <- ifelse(df==0,NA,1-pchisq(Q,df))
    res <- list(data.name=as.character(scoretest),
                statistic = Q, parameter = df,
                p.value=pQ, method = "Score test")
    class(res) <- "htest"
    return(res)    
  }

  ### Likelihood ratio test
  objects <- list(object,...)
  if (length(objects)<2) {
    L0 <- logLik(object)
    L1 <- satmodel(object,logLik=TRUE)
    df <- attributes(L1)$df-attributes(L0)$df; names(df) <- "df"
    ##    Q <- -2*(L0-L1);
    Q <- abs(2*(L0-L1));
    attributes(Q) <- NULL; names(Q) <- "chisq";
    pQ <- ifelse(df==0,NA,1-pchisq(Q,df))

    values <- c(L0,L1); names(values) <- c("log likelihood (model)", "log likelihood (saturated model)")
    res <- list(statistic = Q, parameter = df,
                p.value=pQ, method = "Likelihood ratio test",
                estimate = values)
    class(res) <- "htest"
    return(res)    
  }
  if (length(objects)==2)
    return(comparepair(objects[[1]],objects[[2]]))  
  res <- list()
  for (i in 1:(length(objects)-1)) {
    res <- c(res, list(comparepair(objects[[i]],objects[[i+1]])))
  }
    return(res)
}
