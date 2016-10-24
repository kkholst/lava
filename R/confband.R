##' Add Confidence limits bar to plot
##'
##' @title Add Confidence limits bar to plot
##' @param x Position (x-coordinate if vert=TRUE, y-coordinate otherwise)
##' @param lower Lower limit (if NULL no limits is added, and only the
##' center is drawn (if not NULL))
##' @param upper Upper limit
##' @param center Center point
##' @param line If FALSE do not add line between upper and lower bound
##' @param delta Length of limit bars
##' @param centermark Length of center bar
##' @param pch Center symbol (if missing a line is drawn)
##' @param blank If TRUE a white ball is plotted before the center is
##' added to the plot
##' @param vert If TRUE a vertical bar is plotted. Otherwise a horizontal
##' bar is used
##' @param polygon If TRUE polygons are added between 'lower' and 'upper'.
##' @param step Type of polygon (step-function or piecewise linear)
##' @param ... Additional low level arguments (e.g. col, lwd, lty,...)
##' @seealso \code{confband}
##' @export
##' @keywords iplot
##' @aliases confband forestplot
##' @author Klaus K. Holst
##' @examples
##' plot(0,0,type="n",xlab="",ylab="")
##' confband(0.5,-0.5,0.5,0,col="darkblue")
##' confband(0.8,-0.5,0.5,0,col="darkred",vert=FALSE,pch=1,cex=1.5)
##' 
##' set.seed(1)
##' K <- 20
##' est <- rnorm(K)
##' se <- runif(K,0.2,0.4)
##' x <- cbind(est,est-2*se,est+2*se,runif(K,0.5,2))
##' x[c(3:4,10:12),] <- NA
##' rownames(x) <- unlist(lapply(letters[seq(K)],function(x) paste(rep(x,4),collapse="")))
##' rownames(x)[which(is.na(est))] <- ""
##' signif <- sign(x[,2])==sign(x[,3])
##' forestplot(x,text.right=FALSE)
##' forestplot(x[,-4],sep=c(2,15),col=signif+1,box1=TRUE,delta=0.2,pch=16,cex=1.5)
##' forestplot(x,vert=TRUE,text=FALSE)
##' 
##' z <- seq(10)
##' zu <- c(z[-1],10)
##' plot(z,type="n")
##' confband(z,zu,rep(0,length(z)),col=Col("darkblue"),polygon=TRUE,step=TRUE)
##' confband(z,zu,zu-2,col=Col("darkred"),polygon=TRUE,step=TRUE)
##' 
##' z <- seq(0,1,length.out=100)
##' plot(z,z,type="n")
##' confband(z,z,z^2,polygon="TRUE",col=Col("darkblue"))
##' 
##' set.seed(1)
##' k <- 10
##' x <- seq(k)
##' est <- rnorm(k)
##' sd <- runif(k)
##' val <- cbind(x,est,est-sd,est+sd)
##' par(mfrow=c(1,2))
##' plot(0,type="n",xlim=c(0,k+1),ylim=range(val[,-1]),axes=FALSE,xlab="",ylab="")
##' axis(2)
##' confband(val[,1],val[,3],val[,4],val[,2],pch=16,cex=2)
##' plot(0,type="n",ylim=c(0,k+1),xlim=range(val[,-1]),axes=FALSE,xlab="",ylab="")
##' axis(1)
##' confband(val[,1],val[,3],val[,4],val[,2],pch=16,cex=2,vert=FALSE)
confband <- function(x,lower,upper,center=NULL,line=TRUE,delta=0.07,centermark=0.03,
                     pch,blank=TRUE,vert=TRUE,polygon=FALSE,step=FALSE,...) {
    if (polygon) {
        if (step) {
            x1 <- rep(x,each=2)[-1]
            y1 <- rep(lower, each=2);  y1 <- y1[-length(y1)]
            x2 <- rep(rev(x),each=2); x2 <- x2[-length(x2)]
            y2 <- rep(rev(upper),each=2)[-1]
            xx <- c(x1,x2)
            yy <- c(y1,y2)
        } else {
            xx <- c(x,rev(x))
            yy <- c(lower,rev(upper))
        }
        polygon(xx,yy,...)
        return(invisible(NULL))
    }
    if (vert) {
        ## lower <- lower[length(x)]
        ## upper <- upper[length(x)]
        ## center <- center[length(x)]
        if (line && !missing(lower) && !missing(upper))
            segments(x,lower,x,upper,...)
        if (!missing(lower))
            segments(x-delta,lower,x+delta,lower,...)
        if (!missing(upper))
            segments(x-delta,upper,x+delta,upper,...)
        if (!is.null(center)) {
            if (!missing(pch)) {
                if (blank)
                    points(x,center,pch=16,col="white")
                points(x,center,pch=pch,...)
            } else {
                segments(x-centermark,center,x+centermark,center,...)
            }
        }
    } else {
        if (line && !missing(lower) && !missing(upper))
            segments(lower,x,upper,x,...)
        if (!missing(lower))
            segments(lower,x-delta,lower,x+delta,...)
        if (!missing(upper))
            segments(upper,x-delta,upper,x+delta,...)
        
        if (!is.null(center)) {
            if (!missing(pch)) {
                if (blank)
                    points(center,x,pch=16,col="white")
                points(center,x,pch=pch,...)
            } else {
                ##segments(center,x-centermark,center,x+centermark,...)
            }
        }
    }
    if (missing(lower)) lower <- NULL
    if (missing(upper)) upper <- NULL
    invisible(c(x,lower,upper,center))
}


##' @export
forestplot <- function(x,lower,upper,line=0,labels,text=TRUE,text.right=text,delta=0,axes=TRUE,cex=1,pch=15,xlab="",ylab="",sep,air,xlim,ylim,mar,box1=FALSE,box2=FALSE,vert=FALSE,cex.axis=1,cex.estimate=0.6,...) {
    if (is.matrix(x)) {
        lower <- x[,2]; upper <- x[,3]
        if (ncol(x)>3) cex <- x[,4]
        x <- x[,1]
    }
    if (missing(mar)) {
        if (vert) {            
            mar <- c(8,4,1,1)
        } else {
            mar <- c(4,8,1,1)
        }
    }
    if (missing(labels)) labels <- names(x)
    K <- length(x)
    def.par <- par(no.readonly=TRUE)
    on.exit(par(def.par))
    if (text.right) {
        if (vert) {
            layout(rbind(2,1),heights=c(0.2,0.8))
        } else {
            layout(cbind(1,2),widths=c(0.8,0.2))
        }
    } else {
        layout(1)
    }
    if (vert) {
        if (missing(ylim)) {
            if (missing(air)) air <- max(upper-lower,na.rm=TRUE)*0.4
            ylim <- range(c(x,lower-air,upper+air),na.rm=TRUE)
        }
        if (missing(xlim)) xlim <- c(1,K)        
    } else {
        if (missing(ylim)) ylim <- c(1,K)
        if (missing(xlim)) {
            if (missing(air)) air <- max(upper-lower,na.rm=TRUE)*0.4
            xlim <- range(c(x,lower-air,upper+air),na.rm=TRUE)
        }
    }
    par(mar=mar) ## bottom,left,top,right
    plot(0,type="n",axes=FALSE,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...)
    if (box1) box()
    if (axes) {
        if (vert) {
            axis(2,cex.axis=cex.axis)
        } else {
            axis(1,cex.axis=cex.axis)
        }
    }
    if (!is.null(line)) {
        if (vert) {
            abline(h=line,lty=2,col="lightgray")
        } else {
            abline(v=line,lty=2,col="lightgray")
        }
    }
    if (!missing(sep)) {
        if (vert) {
            abline(v=sep+.5,col="gray")
        } else {
            abline(h=sep+.5,col="gray")
        }
    }
    confband(seq(K),lower,upper,x,pch=pch,cex=cex,vert=vert,blank=FALSE,...)
    if (vert) {
        mtext(labels,1,at=seq(K),las=2,line=1,cex=cex.axis)
    } else {
        mtext(labels,2,at=seq(K),las=2,line=1,cex=cex.axis)
    }
    if (text) {
        xpos <- upper
        if (text.right) {
            if (vert) {
                par(mar=c(0,mar[2],0,mar[4]))
            } else {
                par(mar=c(mar[1],0,mar[3],0))
            }
            plot.new()
            if (vert) {
                plot.window(xlim=xlim,ylim=c(0,0.5))
            } else {
                plot.window(ylim=ylim,xlim=c(0,0.5))
            }
            if (box2) box()
            xpos[] <- 0
        }
        for (i in seq_len(K)) {
            st <- paste0(formatC(x[i])," (",formatC(lower[i]),";",formatC(upper[i]),")")
            if (vert) {
                st <- paste(" ", st)
                if (!is.na(x[i])) graphics::text(i,xpos[i],st,xpd=TRUE,las=2,srt=90, offset=0, pos=4, cex=cex.estimate)
            } else {
                if (!is.na(x[i])) graphics::text(xpos[i],i,st,xpd=TRUE,pos=4,cex=cex.estimate)
            }
        }
    }
}
