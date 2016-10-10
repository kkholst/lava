
backdoor <- function(object,f,cond=NULL,...) {   
   y <- getoutcome(f,sep='|')
   x <- attr(y,'x')
   if (length(x)>1) {
     cond <- all.vars(x[[2]])
     x <- all.vars(x[[1]])
   }
   nod <- vars(object) 
   des <- descendants(object,x)   

   ## We only look at back-door paths:   
   ch <- children(object,x)
   g0 <- cancel(object,toformula(x,ch))
   ## Possible set of variables to condition on:
   cset <- setdiff(nod,c(des,x,y))
   if (is.null(cond)) {
       ## Here we could make traverse through back-doors to find possible sets to condition on
   }
   dsep(g0,c(y,x),cond=cond) && !any(cond%in%des)
}
