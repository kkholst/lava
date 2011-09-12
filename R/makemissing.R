
makemissing <- function(data,p=0.2,cols=1:ncol(data),rowwise=FALSE,nafun=function(x) x) {
  p <- rep(p,length.out=length(cols))
  if (!rowwise)
  for (i in 1:length(cols)) {
    data[rbinom(nrow(data),1,p[i])==1,cols[i]] <- NA
  }
  else
    data[which(rbinom(nrow(data),1,p)==1),cols] <- NA
  return(nafun(data))
}
