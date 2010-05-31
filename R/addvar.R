`addvar` <-
function(x,...) UseMethod("addvar")
`addvar<-` <-
function(x,...,value) UseMethod("addvar<-")

`addvar<-.lvm` <-function(x,...,value) {
  if (class(value)[1]=="formula") {
    return(addvar(x,all.vars(value),...))
  }
  addvar(x, var=value, ...)
}

`addvar.lvm` <-
function(x, var, debug=FALSE, silent=FALSE,...) {
  new <- setdiff(var,vars(x))
  Debug(new, debug)
  if (length(new)>0)
    for (i in new) {
      Graph(x) <- addNode(i, Graph(x))
      x <- addattr(x,attr="shape",var=i,val="rectangle")
      N <- nrow(x$cov)
      if (is.null(N)) {
        N <- 0
        x$cov <- matrix(1); x$covfix <- x$fix <- x$par <- x$covpar <- matrix(NA)
        x$mean <- list(NA)
      } else {
        x$par <- rbind(cbind(x$par, rep(NA,N)), rep(NA,N+1)); ## Add regression labels
        x$covpar <- rbind(cbind(x$covpar, rep(NA,N)), rep(NA,N+1)); ## Add covariance labels
        x$cov <- rbind(cbind(x$cov, rep(0,N)), rep(0,N+1)); ## Add covariance
        x$fix <- rbind(cbind(x$fix, rep(NA,N)), rep(NA,N+1)); ##
        x$covfix <- rbind(cbind(x$covfix, rep(NA,N)), rep(NA,N+1)); ##
        x$mean <- c(x$mean, NA)
      }
      x$cov[N+1,N+1] <- 1
      names(x$mean)[N+1] <- 
        colnames(x$covfix)[N+1] <- rownames(x$covfix)[N+1] <-
          colnames(x$fix)[N+1] <- rownames(x$fix)[N+1] <-
            colnames(x$covpar)[N+1] <- rownames(x$covpar)[N+1] <-               
              colnames(x$par)[N+1] <- rownames(x$par)[N+1] <- 
                colnames(x$cov)[N+1] <- rownames(x$cov)[N+1] <- i
      myexpr <- paste("c(",i,"=expression(",i,"))", sep="\"")
      labels(x) <- (eval(parse(text=myexpr)))
      if (!silent)
        cat("\tAdded '", i, "' to model.\n", sep="")
    }
  return(x)
}
