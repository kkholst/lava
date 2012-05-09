##' @export
`nodecolor<-` <-
function(object,var,...,value) UseMethod("nodecolor<-")

##' @S3method nodecolor<- default
`nodecolor<-.default` <-
  function(object, var=vars(object), border, labcol, shape, lwd, ..., value) {
    if (length(var)>0 & length(value)>0) {
      if (class(var)[1]=="formula") var <- all.vars(var)
      object <- addattr(object,attr="fill",var=var,val=value)
      if (!missing(border))
        object <- addattr(object,attr="col",var=var,val=border)
      if (!missing(shape))
        object <- addattr(object,attr="shape",var=var,val=shape)
      if (!missing(labcol))
        object <- addattr(object,attr="textCol",var=var,val=labcol)
      if (!missing(lwd))
        object <- addattr(object,attr="lwd",var=var,val=lwd)

    }
    return(object)
  }


