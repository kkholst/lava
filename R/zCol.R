##' This function transforms a standard color (e.g. "red") into an
##' transparent RGB-color (i.e. alpha-blend<1). 
##'
##' This only works for certain graphics devices (Cairo-X11 (x11 as of R>=2.7), quartz, pdf, ...).
##' @title Generate a transparent RGB color
##' @param col Color (numeric or character)
##' @param alpha Degree of transparency (0,1)
##' @param locate Choose colour (with mouse) 
##' @return   A character vector with elements of 7 or 9 characters, '"\#"'
##'  followed by the red, blue, green and optionally alpha values in
##' hexadecimal (after rescaling to '0 ... 255').
##' @author Klaus K. Holst
##' @examples
##' plot(runif(1000),cex=runif(1000,0,4),col=Col(c("darkblue","orange"),0.5),pch=16)
##' @keywords color
##' @export
Col <- function(col,alpha=0.2,locate=0) {
    if (locate>0) {
        ytop    <- rep(seq(1/26,1,by=1/26),each=26)[1:657]
        ybottom <- rep(seq(0,1-1/26,by=1/26),each=26)[1:657]
        xleft   <- rep(seq(0,1-1/26,by=1/26),times=26)[1:657]
        xright  <- rep(seq(1/26,1,by=1/26),times=26)[1:657]
        pall    <- round(col2rgb(colors())/256)
        pall    <- colSums(pall) ; pall2 <- character(0)
        pall2[pall>0]   <- "black"
        pall2[pall==0]  <- "white"
        
        par(mar=c(0,0,1,0))
        
        plot.new()
        title(main="Palette of colors()")
        rect(xleft,ybottom,xright,ytop,col=colors())
        text(x=xleft+((1/26)/2)
             ,y=ytop-((1/26)/2)
             ,labels = 1:657
             ,cex=0.55
             ,col=pall2)
        
        
        colmat    <- matrix(c(1:657,rep(NA,26^2-657)),byrow=T,ncol=26,nrow=26)
        cols        <- NA
        i        <- NA
        for(i in 1:locate)
            {
                h    <- locator(1)
                if(any(h$x<0,h$y<0,h$x>1,h$y>1)) stop("locator out of bounds!")
                else {
                    cc        <- floor(h$x/(1/26))+1
                    rr        <- floor(h$y/(1/26))+1            
                    cols[i]    <- colors()[colmat[rr,cc]]
                } 
            } 
        return(cols)
    } 
    
    sapply(col, function(x)
           do.call(rgb,as.list(c(col2rgb(x)/255,alpha)))
           )
}
