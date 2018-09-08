require(mets)
require(lava)

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
