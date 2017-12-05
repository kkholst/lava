toTheta <- function(mu,Sigma,p) {
  theta <- c(as.vector(t(mu)), as.vector(t(Sigma)), p[-nrow(mu)])
  return(theta)
}
toPar <- function(theta, D, k) {
  mus <- Sigmas <- c()
  for (j in 1:k) {
    muj.idx <- (1+((j-1)*D)):(j*D)
    mus <- rbind(mus, theta[muj.idx])
    Sigmaj.start <- k*D+1 + ((j-1)*D^2)
    Sigmaj.idx <- Sigmaj.start + 1:D^2-1
    Sigmas <- rbind(Sigmas, theta[Sigmaj.idx])
  };  ps <- tail(theta,k-1); ps <- c(ps,1-sum(ps))
  return(list(mu=mus, Sigma=Sigmas, p=ps))
}
getMeanVar <- function(object,k,iter,...) {
  if (missing(iter))
    pp <- with(object,toPar(pars,D,k))
  else
    pp <- with(object,toPar(thetas[iter,],D,k))
  res <- list()
  for (i in 1:object$k) {
    mu <- pp$mu[i,]
    V <- matrix(pp$Sigma[i,],ncol=object$D);
    res <- c(res, list(list(mean=mu, var=V)))
  }
  if (missing(k))
    return(res)
  else
    return(res[[k]])  
}



#' Estimate mixture latent variable model
#' 
#' Estimate mixture latent variable model
#' 
#' Estimate parameters in a mixture of latent variable models via the EM
#' algorithm.
#' 
#' @param data \code{data.frame}
#' @param k Number of mixture components
#' @param theta Optional starting values
#' @param steps Maximum number of iterations
#' @param tol Convergence tolerance of EM algorithm
#' @param lambda Added to diagonal of covariance matrix (to avoid
#' singularities)
#' @param mu Initial centres (if unspecified random centres will be chosen)
#' @param silent Turn on/off output messages
#' @param extra Extra debug information
#' @param ... Additional arguments parsed to lower-level functions
#' @return A \code{mixture} object
#' @author Klaus K. Holst
#' @seealso \code{mixture}
#' @keywords models, regression
#' @examples
#' 
#' data(faithful)
#' set.seed(1)
#' M1 <- mvnmix(faithful[,"waiting",drop=FALSE],k=2)
#' M2 <- mvnmix(faithful,k=2)
#' if (interactive()) {
#'     par(mfrow=c(2,1))
#'     plot(M1,col=c("orange","blue"),ylim=c(0,0.05))
#'     plot(M2,col=c("orange","blue"))
#' }
#' 
#' @export mvnmix
mvnmix <- function(data, k=2, theta, steps=500,
                 tol=1e-16, lambda=0,
                 mu=NULL,
                 silent=TRUE, extra=FALSE, ...
                 )  {

  if (k<2) stop("Only one cluster")
  ## theta = (mu1, ..., muk, Sigma1, ..., Sigmak, p1, ..., p[k-1])
  if (is.vector(data)) data <- matrix(data,ncol=1)
  if (is.data.frame(data)) data <- as.matrix(data)  
  i <- 0
  E <- tol
  n <- nrow(data)
  D <- ncol(data)
  yunique <- unique(data)

  if (missing(theta)) {
    mus <- c()
    if (!is.null(mu)) {
      mus <- mu
    } else
    for (j in 1:k) {
      mus <- c(mus, yunique[sample(nrow(yunique),1),])
    }
    Sigmas <- rep(as.vector(cov(data)),k)
    ps <- rep(1/k,k-1)
    theta <- c(mus,Sigmas,ps)
  }

  theta0 <- theta
  if (!silent)
    cat(i,":\t", paste(formatC(theta0),collapse=" "),"\n")
  thetas <- members <- c()
  while ((i<steps) & (E>=tol)) {
    if (extra)
      thetas <- rbind(thetas, theta)
    pp <- toPar(theta,D,k)
    mus <- pp$mu; Sigmas <- pp$Sigma; ps <- pp$p
    ## E(expectation step)
    phis <- c()
    for (j in 1:k) {
        C <- matrix(Sigmas[j,],ncol=D); diag(C) <- diag(C)+lambda ## Assure C is not singular
        phis <- cbind(phis, lava::dmvn(data,mus[j,],C))
    }
    gammas <- c()
    denom <- t(ps%*%t(phis))
    for (j in 1:k) {
      gammas <- cbind(gammas, ps[j]*phis[,j]/denom)
    }
    
    sqrtgammas <- sqrt(gammas)
    ## M(aximization step)
    mus.new <- c()
    Sigmas.new <- c()
    for (j in 1:k) {
      if (!is.null(mu)) mus.new <- mu
      else {
        mu.new <- colSums(gammas[,j]*data)/sum(gammas[,j])
        mus.new <- rbind(mus.new, mu.new)
      }
      ##      browser()
      ##tcrossprod(t(y)-mus.new[[j]])
      wcy <- sqrtgammas[,j]*t(t(data)-mus.new[j,])
      Sigma.new <- t(wcy)%*%wcy/sum(gammas[,j])      
      ## Sigma.new <- 0
      ## for (l in 1:n) {
      ##    Sigma.new <- Sigma.new + gammas[l,j]*(y[l,]-mu.new)%*%t(y[l,]-mu.new)
      ##  }; Sigma.new <- Sigma.new/sum(gammas[,j])      
      Sigmas.new <- rbind(Sigmas.new, as.vector(Sigma.new))
    }; ps.new <- colMeans(gammas)
    theta.old <- theta
    if (extra)
      members <- cbind(members,
                       apply(gammas,1,function(x) order(x,decreasing=TRUE)[1]))
  theta <- toTheta(mus.new,Sigmas.new,ps.new)
  E <- sum((theta-theta.old)^2)
    i <- i+1
    iter <- i    
    if (!silent)
      cat(i,":\t", paste(formatC(theta),collapse=" "),
          ",\t\te=",formatC(E), "\n",sep="")
  }

  myvars <- colnames(data)
  if (is.null(myvars)) myvars <- colnames(data) <- paste("y",1:NCOL(data),sep="")
  data <- as.data.frame(data)
  m <- lvm(myvars,silent=TRUE); m <- covariance(m,myvars,pairwise=TRUE)
  models <- datas <- c()
  for (i in 1:k) {
    models <- c(models, list(m))
    datas <- c(datas, list(data))
  }

  membership <- apply(gammas,1,function(x) order(x,decreasing=TRUE)[1])
  res <- list(pars=theta, thetas=thetas , gammas=gammas, member=membership,
              members=members, k=k, D=D, data=data, E=E,
              prob=rbind(colMeans(gammas)),
              iter=iter,
              models=models,      
              multigroup=multigroup(models,datas)              
              )
  class(res) <- c("mvn.mixture","lvm.mixture")

  parpos <- c()
  npar1 <- D+D*(D-1)/2  
  for (i in 1:k)
    parpos <- c(parpos, list(c(seq_len(D)+(i-1)*D, k*D + seq_len(npar1)+
                               (i-1)*(npar1))))
  
  theta <- c(unlist(lapply(getMeanVar(res),function(x) x$mean)),
             unlist(lapply(getMeanVar(res),function(x) c(diag(x$var),unlist(x$var[upper.tri(x$var)])))))
  res$theta <- rbind(theta)
  res$parpos <- parpos
  res$opt <- list(estimate=theta)
  res$vcov <- solve(information(res,type="E"))   
  return(res)
}


##' @export
print.mvn.mixture <- function(x,...) {
  par <- toPar(x$pars,x$D,x$k)
  space <- paste(rep(" ",12),collapse="")
  for (i in 1:x$k) {
    cat("Cluster ",i," (p=",formatC(par$p[i]),"):\n",sep="")
    cat(rep("-",50),"\n",sep="")
    cat("\tcenter = \n ",space,paste(formatC(par$mu[i,]),collapse=" "),sep="")
    cat("\n\tvariance = \t");
    V <- matrix(formatC(par$Sigma[i,],flag=" "),ncol=x$D);
    colnames(V) <- rep("",x$D); rownames(V) <- rep(space,x$D)
    print(V, quote=FALSE)
    cat("\n")   
  }
  invisible(par)
}

##' @export
plot.mvn.mixture <- function(x, label=2,iter,col,alpha=0.5,nonpar=TRUE,...) {
  opts <- list(...)
  ##  cols <- opts$col; if(is.null(cols)) cols <- 1:gmfit$k
  if (missing(col)) col <- 1:x$k
  lwd <- opts$lwd; if (is.null(lwd)) lwd <- 2
  cex <- opts$cex; if(is.null(cex)) cex <- 0.9
  y <- as.matrix(x$data)
  if (is.vector(y)) y <- matrix(y,ncol=1)
  pp <- getMeanVar(x,iter=iter)
  D <- ncol(y)
  pi <- colSums(x$gammas)/nrow(x$gammas)

  if (D==1) {
    if (nonpar)      
      plot(density(as.vector(y)), main="", ...)
    else
      plot(density(y), main="", type="n", col="lightgray", ...)
    if (!is.null(label)) {
      for (i in 1:x$k) {
        rug(y[x$member==i], col=col[i])
      }
    }
    else
      rug(y)
    cc <- par("usr")
    {
      mycurve <- function(xx) {
        a <- 0;
        for (i in 1:(x$k)) 
          a <- a+pi[i]*dnorm(xx,pp[[i]]$mean,sqrt(pp[[i]]$var[1]))
        a
      }
      curve(mycurve, from=cc[1], to=cc[2], add=TRUE, lwd=lwd,...)
    }
  }
  if (D==2) {
      if (!requireNamespace("ellipse")) stop("ellipse required")
    plot(y, type="n", ...)

    for (i in 1:x$k) {
      C1 <- with(pp[[i]], ellipse::ellipse(var, centre=mean))
      lines(C1, col=col[i], lwd=lwd)      
    }
    
    if (!is.null(label)) {
      for (i in 1:x$k) {
        if (label==1 | missing(iter)) {
          pot <- y[which(x$member==i),]
        }
        else {
          pot <- y[which(x$members[,iter]==i),]
        }
        points(pot, cex=cex, pch=16, col=do.call(rgb, as.list(col2rgb(col[i])/255,alpha)))
      }
    }
    else
      points(y, cex=cex)
  }
  if (D==3) {    
    if (!requireNamespace("rgl")) stop("rgl required")
    rgl::plot3d(y, type="n", box=FALSE)
    for (i in 1:x$k) {
        pot <- y[which(x$member==i),]
        rgl::plot3d(pot, type="s", radius=0.1, col=col[i], add=TRUE)
        ee <- rgl::ellipse3d(pp[[i]]$var,centre=pp[[i]]$mean)
        rgl::plot3d(ee, col=col[i], alpha=alpha, add = TRUE)      
    }
  }  
}
