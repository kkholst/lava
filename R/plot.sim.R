##' @export
##' @export density.sim
density.sim <- function(x,...,plot.type="single") {
    plot.sim(x,...,scatter.plot=FALSE,plot.type=plot.type)
}

##' Plot method for simulation 'sim' objects
##'
##' Density and scatter plots
##' @examples
##' n <- 1000
##' val <- cbind(est1=rnorm(n,sd=1),est2=rnorm(n,sd=0.2),est3=rnorm(n,1,sd=0.5),
##'              sd1=runif(n,0.8,1.2),sd2=runif(n,0.1,0.3),sd3=runif(n,0.25,0.75))
##' 
##' plot.sim(val,estimate=c(1,2),true=c(0,0),se=c(4,5),equal=TRUE,scatter.plot=TRUE)
##' plot.sim(val,estimate=c(1,3),true=c(0,1),se=c(4,6),xlim=c(-3,3),
##' 	scatter.ylim=c(-3,3),scatter.plot=TRUE)
##' plot.sim(val,estimate=c(1,2),true=c(0,0),se=c(4,5),equal=TRUE,
##' 	plot.type="single",scatter.plot=TRUE)
##' plot.sim(val,estimate=c(1),se=c(4,5,6),plot.type="single",scatter.plot=TRUE)
##' plot.sim(val,estimate=c(1,2,3),equal=TRUE,scatter.plot=TRUE)
##' plot.sim(val,estimate=c(1,2,3),equal=TRUE,byrow=TRUE,scatter.plot=TRUE)
##' plot.sim(val,estimate=c(1,2,3),plot.type="single",scatter.plot=TRUE)
##' plot.sim(val,estimate=1,se=c(3,4,5),plot.type="single",scatter.plot=TRUE)
##' 
##' density.sim(val,estimate=c(1,2,3),density=c(0,10,10),angle=c(0,45,-45))
##' @aliases density.sim plot.sim
##' @export
##' @export plot.sim
##' @param x sim object
##' @param estimate columns with estimates
##' @param se columns with standard error estimates
##' @param true (optional) vector of true parameter values
##' @param names (optional) names of estimates
##' @param auto.layout Auto layout (default TRUE)
##' @param byrow Add new plots to layout by row
##' @param type plot type
##' @param ask if TRUE user is asked for input, before a new figure is drawn
##' @param col colour (for each estimate)
##' @param pch plot symbol
##' @param cex point size
##' @param lty line type
##' @param lwd line width
##' @param legend legend 
##' @param legendpos legend position
##' @param cex.legend size of legend text
##' @param plot.type 'single' or 'multiple' (default)
##' @param polygon if TRUE fill the density estimates with colour
##' @param density if non-zero add shading lines to polygon
##' @param angle shading lines angle of polygon
##' @param cex.axis Font size on axis
##' @param alpha Semi-transparent level (1: non-transparent, 0: full)
##' @param main Main title
##' @param cex.main Size of title font
##' @param equal Same x-axis and y-axis for all plots
##' @param delta Controls the amount of space around axis limits
##' @param ylim y-axis limits
##' @param xlim x-axis limits
##' @param ylab y axis label
##' @param xlab x axis label
##' @param rug if TRUE add rug representation of data to x-axis
##' @param rug.alpha rug semi-transparency level
##' @param line.col line colour (running mean, only for scatter plots)
##' @param line.lwd line width (running mean, only for scatter plots)
##' @param line.lty line type (running mean, only for scatter plots)
##' @param line.alpha line transparency 
##' @param scatter.ylab y label for density plots
##' @param scatter.ylim y-axis limits for density plots
##' @param scatter.xlim x-axis limits for density plots
##' @param scatter.alpha semi-transparency of scatter plot
##' @param scatter.col scatter plot colour
##' @param border border colour of density estimates
##' @param true.lty true parameter estimate line type
##' @param true.col true parameter colour
##' @param true.lwd true parameter line width
##' @param density.plot if TRUE add density plot
##' @param scatter.plot if TRUE add scatter plot 
##' @param running.mean if TRUE add running average estimate to scatter plot
##' @param ... additional arguments to lower level functions
plot.sim <- function(x,estimate,se=NULL,true=NULL,
             names=NULL,
             auto.layout=TRUE,
             byrow=FALSE,
             type="p",
             ask=grDevices::dev.interactive(),             
             col=c("gray60","orange","darkblue","seagreen","darkred"),
             pch=16,cex=0.5,lty=1,lwd=0.3,
             legend,
             legendpos="topleft",
             cex.legend=0.8,
             plot.type=c("multiple","single"),
             polygon=TRUE,
             density=0,
             angle=-45,
             cex.axis=0.8,
             alpha=0.2,
             main,
             cex.main=1,
             equal=FALSE,
             delta=1.15,
             ylim=NULL,
             xlim=NULL,
             ylab="",
             xlab="",
             rug=TRUE,
             rug.alpha=0.5,
             line.col=scatter.col,             
             line.lwd=1,
             line.lty=1,
             line.alpha=1,
             scatter.ylab="Estimate",
             scatter.ylim=NULL,
             scatter.xlim=NULL,
             scatter.alpha=0.5,
             scatter.col=col,
             border=col,
             true.lty=2,true.col="gray70",true.lwd=1.2,
             density.plot=TRUE,
             scatter.plot=FALSE,
             running.mean=scatter.plot,
             ...) {

    if (missing(estimate)) {
        estimate <- seq(ncol(x))
    }
    if (is.null(estimate)) {
        av <- apply(x[,drop=FALSE],2,function(z) cumsum(z)/seq(length(z)))
        graphics::matplot(x,type="p",pch=pch, cex=cex, col=col,...)
        graphics::matlines(av,type="l",col=col,lty=lty,...)
        if (!is.null(true)) abline(h=true,lty=true.lty,...)
        if (missing(legend)) legend <- colnames(x)
        if (!is.null(legend))
            graphics::legend(legendpos,legend=legend,bg="white",
                             col=scatter.col,lty=lty,pch=pch,...)
        return(invisible(NULL))
    }
    if (is.character(estimate)) {
        estimate <- match(estimate,colnames(x))
    }

    K <- length(estimate)
    est <- tru <- c()
    if (length(se)>0) {
        if (K==1 && !is.list(se))
            se <- list(se)
        else se <- as.list(se)
    } else {
        est <- estimate; tru <- true
    }
    for (i in seq_along(estimate)) {
        est <- c(est,list(rep(estimate[i],length(se[[i]]))))
        if (!is.null(true)) tru <- c(tru,list(rep(true[i],length(se[[i]]))))
    }

    if (length(se)>0) {
        for (i in seq_along(se)) {
            if (is.character(se[[i]])) se[[i]] <- match(se[[i]],colnames(x))
        }

    }
    ss <- summary.sim(x,estimate=unlist(est),se=unlist(se),true=unlist(tru),names=names)

    oldpar <- NULL
    on.exit({
        par(oldpar)
        return(invisible(ss))
    })

    single <- tolower(plot.type[1])=="single"

    if (auto.layout) {
        nc <- (scatter.plot || running.mean) + density.plot
        nr <- min(6,K)
        if (single) nr <- 1
        oma.multi = c(2, 0, 2, 0)
        mar.multi = c(1.5, 4.1, 1, 1)
        oldpar <- par(mar=mar.multi, oma=oma.multi,
                      cex.axis=cex.axis,las=1,
                      ask=FALSE)
        if (byrow) {
            par(mfrow=c(nr,nc))
        } else {
            par(mfcol=c(nc,nr))
        }
    }

    dys <- c()
    maxdy <- 0
    if (density.plot)
        for (i in seq(K)) {
            ii <- estimate[i]
            y <- as.vector(x[,ii])
            dy <- stats::density(y,...)
            dys <- c(dys,list(dy))
            maxdy <- max(maxdy,dy$y)
        }

    if (equal || single) {
        if (is.null(scatter.ylim)) {
            rg <- range(x[,estimate])
            rg <- rg+c(-1,1)*abs(diff(rg)*(delta-1))
            scatter.ylim <-  rep(list(rg),K)
        }
        if (density.plot) {
            if (is.null(ylim)) ylim <- rep(list(c(0,maxdy*delta)),K)
            if (is.null(xlim)) xlim <- scatter.ylim
        }
    }

    if (!is.null(ylim)) {
        if (!is.list(ylim)) ylim <- list(ylim)
        ylim <- rep(ylim,length.out=K)
    }
    ylab <- rep(ylab,length.out=K)
    if (!is.null(scatter.ylim)) {
        if (!is.list(scatter.ylim)) scatter.ylim <- list(scatter.ylim)
        scatter.ylim <- rep(scatter.ylim,length.out=K)
    }
    if (!is.null(xlim)) {
        if (!is.list(xlim)) xlim <- list(xlim)
        xlim <- rep(xlim,length.out=K)
    }
    if (missing(main)) {
        main <- NULL
        if (!missing(names)) main <- names
        else if (K>1 && !single) main <- colnames(ss)
    }
    if (!is.null(main)) main <- rep(main,length.out=K)
    if (missing(lty)) {
        lty <- rep(1,K)
        if (single || !polygon) {
            lty <- 1:20
        }
    }

    my.scatter.sim <- function(i,add=FALSE,colors,...) {
        ii <- estimate[i]
        if (!missing(colors)) {
            scatter.col <- line.col <- true.col <- colors[1]
        }
        y <- as.vector(x[,ii])
        args <- list(y,ylab=scatter.ylab[i],col=Col(scatter.col[1],scatter.alpha),cex=cex,pch=pch,type=type)
        if (!is.null(scatter.ylim)) args <- c(args,list(ylim=scatter.ylim[[i]]))
        if (scatter.plot) {
            if (!add) {
                do.call(graphics::plot,args)
            } else {
                do.call(graphics::points,args)
            }
        }
        if (running.mean) {
            lines(cumsum(y)/seq_along(y),col=line.col[1],lwd=line.lwd,lty=line.lty)
            if (!is.null(true))
                abline(h=true[i],lty=true.lty,col=true.col[1],lwd=true.lwd)
        }
    }

    my.density.sim <- function(i,add=FALSE,colors,
                               alphas=alpha,
                               auto.legend=TRUE,
                               densities=NULL,
                               angles=angle,
                               ...) {
        ii <- estimate[i]
        y <- as.vector(x[,ii])
        if (!missing(colors)) {
            col <- border <- colors
            col <- true.col <- colors
        }
        if (density.plot) {
            dy <- stats::density(y)
            if (is.null(ylim)) {
                density.ylim0 <- c(0,max(dy$y)*delta)
            } else {
                density.ylim0 <- ylim[[i]]
            }
            if (is.null(xlim)) {
                density.xlim0 <- range(dy$x)
            } else {
                density.xlim0 <- xlim[[i]]
            }
            if (!add) graphics::plot(0,0,type="n",main="",ylab=ylab,xlab=xlab,ylim=density.ylim0,xlim=density.xlim0)
            if (polygon) {
                with(dy, graphics::polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=Col(col[1],alpha=alphas[1]),border=NA,density=densities[1],angle=angles[1]))
                if (!is.null(border)) with(dy, lines(x,y,col=border[1],lty=lty[1],lwd=lwd[1]))
            } else {
                graphics::lines(dy,main="",lty=lty[1],col=col[1],lwd=lwd[1])
            }
            if (rug) graphics::rug(y,col=Col(col[1],rug.alpha[1]))
            if (!is.null(main) && !(running.mean || scatter.plot)) {
                title(main[i],cex.main=cex.main)
            }
            if (!is.null(true)) {
                abline(v=true[i],lty=true.lty,col=true.col,lwd=true.lwd)
            }
            if (!is.null(se)) {
                se.pos <- match(se[[i]],unlist(se))
                ns <- length(se.pos)+1
                se.alpha <- rep(alphas,length.out=ns)[-1]
                se.border <- rep(border,length.out=ns)[-1]
                se.col <- rep(col,length.out=ns)[-1]
                se.lty <- rep(lty,length.out=ns)[-1]
                se.lwd <- rep(lwd,length.out=ns)[-1]
                xx <- dy$x
                for (j in seq_along(se.pos)) {
                    if (polygon) {
                        yy <- dnorm(xx,mean=ss["Mean",se.pos[j]],sd=ss["SE",se.pos[j]])
                        if (se.alpha[j]>0) graphics::polygon(c(xx,rev(xx)),c(yy,rep(0,length(yy))),col=Col(se.col[j],alpha=se.alpha[j]),border=NA,density=densities[j],angle=angles[j])
                        if (!is.null(border)) lines(xx,yy,col=se.border[j],lty=se.lty[j],lwd=se.lwd[j])
                    } else {
                        graphics::curve(dnorm(x,mean=ss["Mean",se.pos[j]],sd=ss["SE",se.pos[j]]),lwd=se.lwd[j],lty=se.lty[j],col=se.col[j],add=TRUE)
                    }
                 }
                if (auto.legend) legend <- c("Kernel",colnames(ss)[se.pos])
                if (!is.null(legend)) {
                    if (polygon) {
                        dcol <- c(col[1],se.col)
                        fill <- Col(dcol,alpha)
                        fill[which(alpha==0)] <- NA
                        border[which(alpha==0)] <- NA
                        dcol[which(alpha!=0)] <- NA
                        graphics::legend(legendpos,legend,
                                         fill=fill,
                                         lty=1,
                                         col=dcol,
                                         border=border,
                                         cex=cex.legend)
                    } else {
                        graphics::legend(legendpos,legend,
                                         col=c(col[1],se.col),
                                         lty=c(lty[1],se.lty),
                                         lwd=c(lwd[1],se.lwd),
                                         cex=cex.legend)
                    }
                }
            }

        }
    }

    if (single) {
        nk <- unlist(lapply(se,length))
        col <- rep(col,length.out=K)
        for (i in seq(K)) {
            my.scatter.sim(i,add=(i>1),colors=col[i])
        }
        if (!is.null(main) && !byrow) {
            title(main[1],cex.main=cex.main)
        }
        if (missing(legend)) legend <- colnames(x)[estimate]
        legendold <- legend
        legend <- NULL
        alpha <- rep(alpha,length.out=K)
        density <- rep(density,length.out=K)
        angle <- rep(angle,length.out=K)
        for (i in seq_len(K)) {
            alphas <- if (K==1L) alpha else alpha[i]
            densities <- density[i]
            if (!is.null(densities) && densities<1) densities <- NULL
            if (length(se)>0) alphas <- c(alphas,rep(0,nk[i]))
            my.density.sim(i,add=(i>1),colors=col[i],alphas=alphas,
                           densities=densities,
                           angles=angle[i],
                           auto.legend=FALSE)
        }
        if (!is.null(legendold)) {
            legend <- rep(legendold,length.out=K)
            graphics::legend(legendpos,legend,
                             fill=Col(col,alpha),border=col,cex=cex.legend)
        }

    } else {
        for (i in seq(K)) {
            my.scatter.sim(i)
            if (!is.null(main) && !byrow && scatter.plot) {
                title(main[i],cex.main=cex.main)
            }
            my.density.sim(i,auto.legend=missing(legend))
            if (i==1 && ask) par(ask=ask)
        }
    }

}

