##' Plot regression line (with interactions) and partial residuals.
##'
##' @title Plot regression lines
##' @param model Model object (e.g. \code{lm})
##' @param var1 predictor (Continuous or factor)
##' @param var2 Factor that interacts with \code{var1}
##' @param data data.frame to use for prediction (model.frame is used as default)
##' @param ci.lty Line type for confidence limits
##' @param ci Boolean indicating wether to draw pointwise 95\% confidence limits
##' @param level Level of confidence limits (default 95\%)
##' @param pch Point type for partial residuals
##' @param lty Line type for estimated regression lines
##' @param lwd Line width for regression lines
##' @param npoints Number of points used to plot curves
##' @param xlim Range of x axis
##' @param col Color (for each level in \code{var2})
##' @param colpt Color of partial residual points
##' @param alpha Alpha level
##' @param cex Point size 
##' @param delta For categorical \code{var1}
##' @param centermark For categorical \code{var1}
##' @param jitter For categorical \code{var1}
##' @param cidiff For categorical \code{var1}
##' @param mean For categorical \code{var1}
##' @param legend Boolean (add legend) 
##' @param trans Transform estimates (e.g. exponential)
##' @param partres Boolean indicating whether to plot partial residuals
##' @param partse .
##' @param labels Optional labels of \code{var2}
##' @param \dots additional arguments to lower level functions
##' @return list with following members:
##' \item{x}{Variable on the x-axis (\code{var1})}
##' \item{y}{Variable on the y-axis (partial residuals)}
##' \item{predict}{Matrix with confidence limits and predicted values}
##' @author Klaus K. Holst
##' @seealso \code{termplot}
##' @export
##' @examples
##' n <- 100
##' x0 <- rnorm(n)
##' x1 <- seq(-3,3, length.out=n)
##' x2 <- factor(rep(c(1,2),each=n/2), labels=c("A","B"))
##' y <- 5 + 2*x0 + 0.5*x1 + -1*(x2=="B")*x1 + 0.5*(x2=="B") + rnorm(n, sd=0.25)
##' dd <- data.frame(y=y, x1=x1, x2=x2)
##' lm0 <- lm(y ~ x0 + x1*x2, dd)
##' plotConf(lm0, var1="x1", var2="x2")
##' ### points(5+0.5*x1 -1*(x2=="B")*x1 + 0.5*(x2=="B") ~ x1, cex=2)
##' 
##' data(iris)
##' l <- lm(Sepal.Length ~ Sepal.Width*Species,iris)
##' plotConf(l,var2="Species")
##' @keywords hplot, regression
##' 
plotConf <- function(model, ##var1=attributes(terms(model), "terms.label")[2],
                     var1=all.vars(formula(model))[2],
                     var2=NULL,
                     data=NULL,
                     ci.lty=0,
                     ci=TRUE, level=0.95,
                     pch=16, lty=1, lwd=2,
                     npoints=100,
                     xlim,
                     col=NULL,
                     colpt,
                     alpha=0.2,
                     cex=1,
                     delta=0.07,
                     centermark=0.03,
                     jitter=0.2,
                     cidiff=FALSE,
                     mean=TRUE,
                     legend=ifelse(is.null(var1),FALSE,"topright"),
                     trans=function(x) {x},
                     partres=TRUE,
                     partse=FALSE,
                     labels,
                     ...) {
  
  intercept <- coef(model)["(Intercept)"]
  if (is.na(intercept) & !is.null(var2)) {
    model <- update(model,.~.+1)
    warning("Refitted model with an intercept term")
    intercept <- coef(model)["(Intercept)"]  
  }

  dots <- list(...)
  

  if (is.null(data)) {
    curdata <- model.frame(model)
  } else {
    curdata <- get_all_vars(formula(model), data)
  }
  responseorig <- colnames(curdata)[1]

  response <- all.vars(formula(model))[1]
  cname <- colnames(curdata)[-1]
  if (!is.factor(curdata[,var2]) & !is.null(var2)) {
    message("var2 should be a factor (or NULL)")
    curdata[,var2] <- as.factor(curdata[,var2])
    colnames(curdata)[1] <- response
    model <- update(model,as.formula(paste(response,"~.")),data=curdata)
##    thelevels <- unique(curdata[,var2])    
  }
  thelevels <- levels(curdata[,var2])
  if (missing(labels)) labels <- thelevels
  
  k <- ifelse(is.null(var2),1,length(thelevels))
  curpal <- palette()
  if (is.null(col)) { ## assign colors
    col <- seq(k)
    mypal()
  }
  
  if (missing(xlim)) {
    if (!is.null(var1))
      xlim <- range(curdata[,var1])
    else
      xlim <- c(0,length(thelevels))+0.5
  }
  dots$xlim <- xlim
  x <- seq(xlim[1], xlim[2], length.out=npoints)  
  xx <- c()
 

  newdata <- data.frame(id=1:npoints)
  partdata <- curdata[,-1,drop=FALSE]
  for (nn in cname) {
    v <- curdata[,nn]
    if (!is.null(var1) && nn==var1) {
      newdata <- cbind(newdata, rep(x, k))
      partdata[,nn] <- 0
    } else {
      if (is.factor(v)) {
        if (nn%in%var2) {
          newdata <- cbind(newdata, factor(rep(levels(v), each=npoints),levels=thelevels))
          partdata[,nn] <- factor(rep(levels(v)[1], nrow(partdata)))
        } else {
          newdata <- cbind(newdata, factor(rep(levels(v)[1], k*npoints),
                                           levels=levels(v)))
##          partdata[,nn] <- factor(rep(levels(v)[1], nrow(paprtdata)))
          
        }
        ##  levels(newdata[,ncol(newdata)]) <- levels(v)
        ###newdata[,ncol(newdata)] <- relevel(newdata[,ncol(newdata)],levels(v))        
      } else {
        ## ref <- thelevels[1]
        ## if (is.numeric(curdata[,var2])) ref <- 0
        ## if (nn%in%var2) {
        ##   newdata <- cbind(newdata, rep(thelevels, each=npoints))
        ##   partdata[,nn] <- rep(ref,nrow(partdata))
        ## } else {
          newdata <- cbind(newdata, rep(0,k*npoints))         
##        }
        ##        partdata[,nn] <- 0
      }
    }
  };
  colnames(newdata) <- c("_id", cname)
  ##deparse(formula(model)[[2]])

  newdata[,response] <- 0
  Y <- model.frame(model)[,1]
  if(class(Y)=="Surv") Y <- Y[,1] 
  
  if (!any(c("lm","glm","survreg")%in%class(model)[1])) {
    names(newdata)[1] <- response
    if ("mer"%in%class(model)) {
      bb <- lme4::fixef(model)
      XX <- model.matrix(formula(terms(model)),data=newdata)
    }
    else {
      bb <- coef(model)
      XX <- model.matrix(formula(model),data=newdata)
    }
    coefnames <- c("(Intercept)",var1)
    if (!is.null(var2)) {
      var2n <- paste(var2,thelevels[-1],sep="")
      coefnames <- c(coefnames,var2n,paste(var1,var2n,sep=":"),
                     paste(var2n,var1,sep=":"))
    }
    bidx <- which(names(bb)%in%coefnames)
    notbidx <- setdiff(seq(length(bb)),bidx)
    if (class(model)[1]=="geeglm")
      SS <- (summary(model)$cov.unscaled)[bidx,bidx,drop=FALSE]
    else
      SS <- vcov(model)[bidx,bidx,drop=FALSE]
    bb0 <- bb; bb0[notbidx] <- 0

    sqrtSS <- with(eigen(SS),vectors%*%diag(sqrt(values))%*%t(vectors))
    myse <- apply(XX[,bidx,drop=FALSE]%*%sqrtSS,1,function(x) sum(x^2))^0.5    
    ##    myse <- diag((XX[,bidx,drop=FALSE]%*%SS)%*%t(XX[,bidx,drop=FALSE]))^0.5
    ci.all <- list(fit=XX%*%bb,se.fit=myse)
    if (!is.na(bb["(Intercept)"])) bidx <- bidx[-1]
    XX0 <- model.matrix(model); XX0[,bidx] <- 0

    Yhat <- XX0%*%bb
    R <- as.vector(Y-Yhat)

    if ("mer"%in%class(model)) {
      uz <- as.matrix(unlist(lme4:::ranef(model))%*%model@Zt)[1,]
      R <- R-uz
    }
  } else {
    bb <- coef(model)
    partdata[,response] <- 0
    R <- Y-predict(model, newdata=partdata)
    ci.all <- predict(model, newdata=newdata, se.fit=TRUE, interval = "confidence", level=level)
  }
  if (ncol(ci.all$fit)==1) {
    z <- qnorm(1-(1-level)/2)
    ci.all$fit <- cbind(ci.all$fit,ci.all$fit-z*ci.all$se.fit,ci.all$fit+z*ci.all$se.fit)
  }

  intercept <- bb["(Intercept)"]
  pr <- trans(intercept + R)
  
  intercept0 <- 0
  if (is.na(intercept)) {
    intercept <- 0
    if (!is.null(var2)) {
      intercept <- coef(model)[paste(var2,thelevels,sep="")][as.numeric(curdata[,var2])]
      intercept0 <- coef(model)[paste(var2,thelevels[1],sep="")]
    }
  }

  ##Y <- curdata[,response]
  if (is.null(dots$ylim)) {
    if (partres) {
      if (cidiff)
        dots$ylim <- range(pr)
      else
        dots$ylim <- range(trans(c(ci.all$fit)),pr)
    }
    else
      dots$ylim <- trans(range(ci.all$fit))
  }
  if (is.null(dots$ylab))
    dots$ylab <- responseorig

## dots$ylab <- as.expression(dots$ylab)
## dots$xlab <- as.expression(dots$xlab)

  if (is.null(var1)) {
    dots$axes=FALSE
    if (is.null(dots$xlab))
      dots$xlab <- ""    
  } else  {
    if (is.null(dots$xlab))
      dots$xlab <- var1    
  }

  plot.list <- c(x=0,y=0,type="n",dots)
  do.call(plot, plot.list)
  if (is.null(var1)) {
    box()
    axis(2)
    axis(1,at=seq(length(thelevels)),labels)
  }
  
  col.trans <- Col(col,alpha)

  Wrap <- function(k,n) { (1:k-1)%%n +1 }
  col.i <- Wrap(k,length(col));  col.k <- col[col.i]; 
  lty.k <- lty[Wrap(k,length(lty))]
  pch.k <- pch[Wrap(k,length(pch))]

  if (!is.null(var1)) {
   
    for (i in 1:k) {
      ci0 <- trans(ci.all$fit[(npoints*(i-1)+1):(i*npoints),])
      y <- ci0[,1]; yu <- ci0[,3]; yl <- ci0[,2]
      lines(y ~ x, col=col.k[i], lwd=lwd, lty=lty.k[i])

      if (ci) {
        lines(yl ~ x, lwd=1, col=col.k[i], lty=ci.lty)
        lines(yu ~ x, lwd=1, col=col.k[i], lty=ci.lty)
        xx <- c(x, rev(x))
        yy <- c(yl, rev(yu))
        polygon(xx,yy, col=col.trans[col.i[i]], lty=0)
      }
    }

  }
##  if (is.factor(partdata[,var1]))
##    partdata[,cname%nin%c(var1,var2)] <- 0
##  pr <- predict(model, newdata=partdata)
##  b <- coef(model); 
##  X <- model.matrix(model)
##  pr <- X%*%b
  
##  points(curdata[,response] ~ curdata[,var1], col=col[as.numeric(curdata[,var2])], pch=16, cex=0.7)
  ii <- as.numeric(curdata[,var2])

  if (is.null(var1)) {
    xx <- curdata[,var2]
    x <- jitter(as.numeric(xx),jitter)
    if (missing(colpt)) colpt <- col[1]
    if (partres>0)
      points(pr ~ x,pch=pch[1], col=colpt[1], ...)
    positions <- seq(k)
    mycoef <- bb[paste(var2,thelevels,sep="")][-1]
    if (class(model)[1]%in%c("lm","glm"))
      myconf <- confint(model)[paste(var2,thelevels,sep="")[-1],,drop=FALSE]    
    else {
      myconf <- matrix(mycoef,ncol=2,nrow=length(mycoef))
      myconf <- myconf + qnorm(0.975)*cbind((diag(as.matrix(SS))[-1])^0.5)%x%cbind(-1,1)
    }
    for (pos in seq(k)) {
      if (cidiff) {
        if (pos>1) {
          ci0 <- trans(intercept+myconf[pos-1,])
          yl <- ci0[1]; yu <- ci0[2]; y <- trans(intercept+mycoef[pos-1])
        } else {
          yu <- yl <- NULL; y <- trans(intercept)
        }
      } else if (partse) {
        y0 <- pr[xx==levels(xx)[pos]]
        ci0 <- confint(lm(y0~1))
        yl <- ci0[1]; yu <- ci0[2]; y <- trans(mean(y0))
      } else {
        ci0 <- trans(ci.all$fit[(npoints*(pos-1)+1):(pos*npoints),])
        y <- ci0[,1]; yu <- ci0[,3]; yl <- ci0[,2]
      }
      if (!mean) y <- NULL
      confband(pos,yl,yu,delta=delta,center=y,centermark=centermark,col=col[1],lwd=lwd[1],lty=lty[1])  
    }
  } else {    
    if (partres) {
      xx <- curdata[,var1]
      if (!missing(colpt)) {
        points(pr ~ xx, col=colpt, cex=cex, pch=pch[1],...)
      } else {
        if (!is.null(var2))
          points(pr ~ xx, col=col.k[ii], pch=pch.k[ii], cex=cex, ...)
        else
          points(pr ~ xx, col=col[1], pch=pch[1], cex=cex,...)
      }
    }
  }
  
  if (k>1 && legend!=FALSE) {
    if (length(lty)>1)
      legend(legend, legend=thelevels, col=col.k, pch=pch.k, bg="white", lty=lty.k,cex=cex)
    else
      legend(legend, legend=thelevels, col=col.k, pch=pch.k, bg="white",cex=cex)
  }

  palette(curpal)  
  invisible(list(x=xx, y=pr, predict=ci.all))
}
