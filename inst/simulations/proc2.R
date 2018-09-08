library(lava)
library(knitr)
library(kableExtra)
options(knitr.table.format = "latex")

idx <- c("Mean","SD","SE","SE/SD","Coverage","RMSE")         
prep <- function(inp,         
         parname="\\(\\bm{\\beta_1=1}\\)",
         digits=3,...) {
    digits <- 3
    a1 <- format(round(t(inp[idx,])*10^digits,0)/10^digits,digits=digits)
    colnames(a1)[which(colnames(a1)=="SE/SD")] <-
        "\\(\\tfrac{\\text{\\bfseries{SE}}}{\\text{\\bfseries{SD}}}\\)"
    a1 <- cbind(rownames(a1), a1)
    rname <- c(rep(parname,nrow(a1)))
    a1 <- cbind(rname, a1)
    colnames(a1)[1] <- ""

    a <- kable(a1, row.names=FALSE, booktabs = TRUE, caption = "My table", escape = FALSE) %>% 
        #add_header_above(c(" ", " ", "colLabel1"=2, "colLabel2"=2)) %>% 
        kable_styling(latex_options = "hold_position") %>%
        row_spec(0, bold=TRUE) %>%
        ##column_spec(0, bold=TRUE) %>%
        collapse_rows(columns = 1)
    res <- gsub("[\\]+cmidrule[\\{][1-9][\\-][1-9][}]","",as.character(a))
    cat(res,"\n")
}

source("proc.R")


##################################################

res <- loadfiles("normal")
cc <- getcoef(res)
gam2 <- structure(cc[[2]],class=c("sim","matrix"))
gam1 <- structure(cc[[1]],class=c("sim","matrix"))

ss2 <- unlist(lapply(res,function(x) mean(x$agq9$score^2)))
ii2 <- which(abs(ss2)<0.5)
est <- c("twostage","twostage.mixture","iv","iv.rob","laplace","agq9")
(res2 <- summary(gam2[ii2,],estimate=est,
                se=paste0("se.",est),
                confint=TRUE,
                true=rep(0.5,length(est))
                ))
(res1 <- summary(gam1[ii2,],estimate=est,
                se=paste0("se.",est),
                confint=TRUE,
                true=rep(1,length(est))
                ))

colnames(res1) <-
    colnames(res2) <-
    c("TS","TS[M2]","2SLS","2SLS[robust]","Laplace","AGQ9")
prep(res1,parname="\\(\\bm{\\beta_1=1}\\)")
prep(res2,parname="\\(\\bm{\\beta_2=0.5}\\)")



res <- loadfiles("mixture2")
cc <- getcoef(res)
gam2 <- structure(cc[[2]],class=c("sim","matrix"))
gam1 <- structure(cc[[1]],class=c("sim","matrix"))

ss2 <- unlist(lapply(res,function(x) mean(x$agq9$score^2)))
ii2 <- which(abs(ss2)<0.5)
ii2 <- seq(500)
est <- c("twostage","twostage.mixture","iv","iv.rob","laplace","agq9")
est <- c("twostage","twostage.mixture","iv.rob")
(res2 <- summary(gam2[ii2,],estimate=est,
                se=paste0("se.",est),
                confint=TRUE,
                true=rep(0.5,length(est))
                ))
(res1 <- summary(gam1[ii2,],estimate=est,
                se=paste0("se.",est),
                confint=TRUE,
                true=rep(1,length(est))
                ))

colnames(res1) <-
    colnames(res2) <-
    c("TS","TS[M2]","2SLS","2SLS[robust]","Laplace","AGQ9")
prep(res1,parname="\\(\\bm{\\beta_1=1}\\)")
prep(res2,parname="\\(\\bm{\\beta_2=0.5}\\)")



##################################################

val <- get(load("data/exp1.rda"))
res <- summary(val, estimate=1:4, se=5:8, true=c(0,.3,0,.3), confint=TRUE)
prep(res[,c(1,3)])
prep(res[,c(2,4)])


val <- get(load("data/uniform1.rda"))
res1 <- summary(val, estimate=(1:7)*2-1, se=(8:14)*2-1, true=rep(1,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
res2 <- summary(val, estimate=(1:7)*2, se=(8:14)*2, true=rep(.5,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
prep(res1)
prep(res2)

val <- get(load("data/uniform1b.rda"))
res1 <- summary(val, estimate=(1:7)*2-1, se=(8:14)*2-1, true=rep(1,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
res2 <- summary(val, estimate=(1:7)*2, se=(8:14)*2, true=rep(.5,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
prep(res1)
prep(res2)


val <- get(load("data/uniform2.rda"))
res1 <- summary(val, estimate=(1:7)*2-1, se=(8:14)*2-1, true=rep(1,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
res2 <- summary(val, estimate=(1:7)*2, se=(8:14)*2, true=rep(.5,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
prep(res1)
prep(res2)

val <- get(load("data/uniform2b.rda"))
res1 <- summary(val, estimate=(1:7)*2-1, se=(8:14)*2-1, true=rep(1,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
res2 <- summary(val, estimate=(1:7)*2, se=(8:14)*2, true=rep(.5,7), confint=TRUE,
                names=c("TS1","TS2","TS3","IV","IVrob","Wall","Wallrob"))
prep(res1)
prep(res2)



