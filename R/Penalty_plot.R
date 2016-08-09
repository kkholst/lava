`plot.plvmfit` <- function(x, lambda = "lambda1.abs", type = NULL,
                               add.line = TRUE, line.size = 2,
                               add.point = TRUE, point.shape = 4, point.size = 2,
                               add.best = TRUE, color.selected = TRUE) {
  
  test.best <- !is.null(x$penalty$lambda1.best)
  
  if(is.null(type)){
    if(test.best){
      type <- "criterion"
    }else{
      type <- "path"
    }
  }
  
  if(type == "path"){
    path <- getPath(x, getCoef = "penalized", getLambda = lambda)
    df.Path <- data.table::melt(path, 
                                measure=names(path)[-1], 
                                value.name = "value", variable.name = "coefficient")
    
    ggPath <- ggplot(df.Path, aes_string(y = "value", 
                                         x = lambda, 
                                         group = "coefficient", 
                                         col = "coefficient")
    )
    if(add.line){ggPath <- ggPath + geom_line(size = line.size)}
    if(add.point){ggPath <- ggPath + geom_point(size = point.size, shape = point.shape)}
  
    if(!is.null(x$penalty$lambda1.best)){
      if(color.selected){
        df.Path$selected <- df.Path$coefficient %in% names(coef(x))
        names.selected <- unique(df.Path$coefficient[df.Path$selected])
        names.Nselected <- unique(df.Path$coefficient[!df.Path$selected])
        n.selected <- length(names.selected)
        n.Nselected <- length(names.Nselected)
        color.selected <- rgb(g = seq(0.3,0.7,length.out = n.selected), 0, 0)
        color.Nselected <- rgb(r = seq(0.3,0.7,length.out = n.Nselected), 0, 0)
        
        color.order <- as.character(unique(df.Path$coefficient))
        color.order[color.order %in% names.selected] <- color.selected
        color.order[color.order %in% names.Nselected] <- color.Nselected
        ggPath <- ggPath + scale_color_manual(values = color.order)
      }
      
      if(add.best){
        ggPath <- ggPath + geom_vline(size = line.size/2, 
                                      xintercept = unlist(getPath(x, names = lambda, row = attr(x$penalty$lambda1.best,"row"))), 
                                      linetype = 2, color = "blue")
        }
      
    }
    
    return(ggPath)
    
  }else if(type == "criterion"){
    
    df <- data.frame(lambda1 = x$penalty$lambda1,
                     criterion = x$penalty$performance,
                     optimum = c("no","yes")[(x$penalty$lambda1 == x$penalty$lambda1.best) + 1])
    names(df)[2] <- attr(x$penalty$performance,"criterion")
    
    ggPerf <- ggplot(df, aes_string(x = "lambda1", y = names(df)[2]))
    if(add.line){ggPerf <- ggPerf + geom_line(size = line.size)}
    if(add.point){ggPerf <- ggPerf + geom_point(size = point.size, aes_string(color = "optimum"))}
    
    return(ggPerf)
  }else{
    stop("type must be \"path\" or \"criterion\" \n")
  }
}
