`baptize` <- function(x,...) UseMethod("baptize")

###{{{ baptize.lvm

baptize.lvm <- function(x,labels,overwrite=FALSE,...) {
  p <- describecoef(x, mean=TRUE)  
  MeanFix <- intfix(x)
  RegFix <- regfix(x)
  CovFix <- covfix(x)
  count <- 0
  for (i in 1:length(p)) {
    p0 <- p[[i]]
    if (attributes(p0)$type=="reg") {
      curfix <- RegFix$values[p0[2],p0[1]]
      curlab <- RegFix$labels[p0[2],p0[1]]
      if (all(is.na(c(curfix,curlab))) | overwrite) {
        count <- count+1
        st <- ifelse(missing(labels),paste("p",count,sep=""),labels[count])
        regfix(x,from=p0[2],to=p0[1]) <- st
      }
    } else if (attributes(p0)$type=="cov") {
      curfix <- CovFix$values[p0[2],p0[1]]
      curlab <- CovFix$labels[p0[2],p0[1]]
      if (all(is.na(c(curfix,curlab))) | overwrite) {
        count <- count+1
        st <- ifelse(missing(labels),paste("p",count,sep=""),labels[count])
##        st <- paste("p",count,sep="")
        covfix(x,p0[2],p0[1],exo=FALSE) <- st
      }
    } else { ## Mean parameter
      curfix <- MeanFix[[p0]]
      if (is.na(curfix) | overwrite) {
        count <- count+1
        st <- ifelse(missing(labels),paste("m",count,sep=""),labels[count])
        intfix(x,p0) <- st        
      }
    }
  }
  return(x)
}

###}}}
