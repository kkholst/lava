library(mets)
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

loadfiles <- function(pattern,idx=NULL) {
    if (is.null(idx)) {
        ff <- list.files("data/", paste0("*_",pattern,".rda"))
        idx <- sort(as.numeric(sapply(ff,function(x) strsplit(x,"_")[[1]][2])))
    }
    res <- lapply(idx,
                  function(x) get(load(paste0("data/res_",x,"_",pattern,".rda"))))
    return(res)   
}

getcoef <- function(res,...) {
    ff <- function(x) {
        if (length(x)>0) return(x) else return(NA)
    }
    pp <- c("twostage","twostage.naive","twostage.naive.rob",
            "twostage.mixture",
            "linear",
            "sas.agq","sas.iv",
           "iv",
           "iv.rob",
           "pmle",
           "wall","wall.rob",
           "laplace","agq9")

    gam1 <- lapply(res,function(x) {
        c(twostage=ff(x$twostage[c("u2~u1_1"),1]),
          twostage.naive=ff(x$twostage.naive[c("u2~u1_1"),1]),
          twostage.naive.rob=ff(x$twostage.naive.robust[c("u2~u1_1"),1]),
          ##twostage.mixture=x$twostage.mixture[c("eta~f1"),1]
          twostage.mixture=if (!is.null(x$twostage.mixture)) 
                               ff(x$twostage.mixture[c("u2~u1_1"),1])
                           else
                               ff(x$twostage[c("u2~u1_1"),1]),
          linear=ff(x$linear$coef["u2~u1",1]),
          sas.agq=ff(x$sas.nsem[15,2]),
          sas.iv=ff(x$sas.iv[2,3]),
          iv=ff(x$iv1[c("y1~f1"),1]),
          iv.rob=ff(x$iv1.rob[c("y1~f1"),1]),
          pmle=ff(x$iv.pmle.robust[c("eta~f1"),1]),
          wall=ff(x$wall[2]),
          wall.rob=ff(x$wall.rob[2]),
          laplace=ff(x$laplace$coef[c("eta2~eta1"),1]),
          agq9=ff(x$agq9$coef[c("eta2~eta1"),1]),
          se.twostage=ff(x$twostage[c("u2~u1_1"),2]),
          se.twostage.naive=ff(x$twostage.naive[c("u2~u1_1"),2]),
          se.twostage.naive.rob=ff(x$twostage.naive.robust[c("u2~u1_1"),2]),
          se.twostage.mixture=if (!is.null(x$twostage.mixture))
                                  ff(x$twostage.mixture[c("u2~u1_1"),2])
                              else
                                  ff(x$twostage[c("u2~u1_1"),2]),
          ##se.twostage.mixture=x$twostage.mixture[c("eta~f1"),2],
          se.linear=ff(x$linear$coef["u2~u1",2]),
          se.sas.agq=ff(x$sas.nsem[15,3]),
          se.sas.iv=ff(x$sas.iv[2,4]),
          se.iv=ff(x$iv1[c("y1~f1"),2]),
          se.iv.rob=ff(x$iv1.rob[c("y1~f1"),2]),
          se.pmle=ff(x$iv.pmle.robust[c("eta~f1"),2]),
          se.wall=NA,
          se.wall.rob=NA,
          se.laplace=ff(x$laplace$coef[c("eta2~eta1"),2]),
          se.agq9=ff(x$agq9$coef[c("eta2~eta1"),2]))
    })
    Gam1 <- matrix(unlist(gam1),byrow=TRUE,nrow=length(gam1))
    colnames(Gam1) <- names(gam1[[1]])
    
    gam2 <- lapply(res,function(x) {
        c(twostage=ff(x$twostage[c("u2~u1_2"),1]),
          twostage.naive=ff(x$twostage.naive[c("u2~u1_2"),1]),
          twostage.naive.rob=ff(x$twostage.naive.robust[c("u2~u1_2"),1]),
          ##twostage.mixture=x$twostage.mixture[c("eta~f1"),1]
          twostage.mixture=if (!is.null(x$twostage.mixture)) 
                               ff(x$twostage.mixture[c("u2~u1_2"),1])
                           else
                               ff(x$twostage[c("u2~u1_2"),1]),
          linear=ff(x$linear$coef["u2~u1",1]),
          sas.agq=ff(x$sas.nsem[16,2]),
          sas.iv=ff(x$sas.iv[3,3]),
          iv=ff(x$iv1[c("y1~f2"),1]),
          iv.rob=ff(x$iv1.rob[c("y1~f2"),1]),
          pmle=ff(x$iv.pmle.robust[c("eta~f2"),1]),
          wall=ff(x$wall[3]),
          wall.rob=ff(x$wall.rob[3]),
          laplace=ff(x$laplace$coef[c("eta2~eta1^2"),1]),
          agq9=ff(x$agq9$coef[c("eta2~eta1^2"),1]),
          se.twostage=ff(x$twostage[c("u2~u1_2"),2]),
          se.twostage.naive=ff(x$twostage.naive[c("u2~u1_2"),2]),
          se.twostage.naive.rob=ff(x$twostage.naive.robust[c("u2~u1_2"),2]),
          se.twostage.mixture=if (!is.null(x$twostage.mixture))
                                  ff(x$twostage.mixture[c("u2~u1_2"),2])
                              else
                                  ff(x$twostage[c("u2~u1_2"),2]),
          ##se.twostage.mixture=x$twostage.mixture[c("eta~f1"),2],
          se.linear=ff(x$linear$coef["u2~u1",2]),
          se.sas.agq=ff(x$sas.nsem[16,3]),
          se.sas.iv=ff(x$sas.iv[3,4]),
          se.iv=ff(x$iv1[c("y1~f2"),2]),
          se.iv.rob=ff(x$iv1.rob[c("y1~f2"),2]),
          se.pmle=ff(x$iv.pmle.robust[c("eta~f2"),2]),
          se.wall=NA,
          se.wall.rob=NA,
          se.laplace=ff(x$laplace$coef[c("eta2~eta1^2"),2]),
          se.agq9=ff(x$agq9$coef[c("eta2~eta1^2"),2]))
        })
    Gam2 <- matrix(unlist(gam2),byrow=TRUE,nrow=length(gam2))
    colnames(Gam2) <- names(gam2[[1]])
    return(list(Gam1,Gam2))
}


##################################################

res <- loadfiles("normalc")
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



res <- loadfiles("mixture2c")
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



