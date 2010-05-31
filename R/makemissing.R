
makemissing <- function(data,p=0.2,cols=1:ncol(data),rowwise=FALSE,nafun=function(x) x) {
  if (!rowwise)
  for (i in cols) {
    data[rbinom(nrow(data),1,p)==1,i] <- NA
  }
  else
    data[rbinom(nrow(data),1,p)==1,cols] <- NA
  return(nafun(data))
}
