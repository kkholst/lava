subset.lvm <- function(x, vars, ...) {
  if (missing(vars))
    return(x)
  if (class(vars)[1]=="formula") vars <- all.vars(vars)
  if (!all(vars%in%vars(x))) stop("Not a subset of model")  
  latentvars <- intersect(vars,latent(x))
  g0 <- subGraph(vars, Graph(x))
  res <- graph2lvm(g0)
  if (length(latentvars)>0)
    latent(res) <- latentvars
  res$cov[vars,vars] <- x$cov[vars,vars]
  ## Fixed parameters:
  res$par[vars,vars] <- x$par[vars,vars]
  res$fix[vars,vars] <- x$fix[vars,vars]
  res$covpar[vars,vars] <- x$covpar[vars,vars]
  res$covfix[vars,vars] <- x$covfix[vars,vars]
  res$mean[vars] <- x$mean[vars]
  index(res) <- reindex(res)
  return(res)
}
