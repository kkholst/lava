comparepair <- function(x1,x2) {
  ##if (class(x1)!="lvmfit" | class(x2)!="lvmfit") stop("'lvmfit' object expected.")
  l1 <- logLik(x1);  l2 <- logLik(x2)
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

compare.default <- function(object,...,par,contrast,null) {
  if (!missing(par)) {
    contrast <- rep(0,length(coef(object)))
    contrast[parpos(Model(object),p=par)] <- 1
    contrast <- diag(contrast)
  }
  if (!missing(contrast)) {
    B <- contrast        
    p <- coef(object)
    if (ncol(B)<length(p)) {
      B <- matrix(0,ncol=length(coef(object)))
      B[parpos(Model(object),p=colnames(contrast))] <- 1      
    }
    if (missing(null)) null <- 0
    Q <- t(B%*%p-null)%*%Inverse(B%*%vcov(object)%*%t(B))%*%(B%*%p-null)
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
    return(res)        
  }
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
