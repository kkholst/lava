library(mets)
library(knitr)
library(kableExtra)
options(knitr.table.format = "latex")

loadall <- function(name="unif200", expr=c(kevala.rmse, ts.rmse, ts0.rmse), idx=1:100,
            names,reduce=TRUE) {
    e <- new.env()
    result <- c()    
    for (i in idx) {
        infile <- paste0("save/res_",name,"_",i,".rda")
        suppressWarnings(re <- try(get(load(infile)),silent=TRUE))
        if (inherits(res,"try-error")) browser()
        if (!inherits(res,"try-error")) {            
            if (is.null(try(eval(substitute(expr)),silent=TRUE))) {
                reduce <- FALSE
                result <- list(res)
            } else {
                for (x in base::names(res)) assign(x,res[[x]],envir=e)
                result <- c(result, list(eval(substitute(expr),envir=e)))
            }
        }
    }
    if (reduce) {
        result <- Reduce(rbind,result)
        if (!missing(names))
            dimnames(result) <- list(NULL, names)
    }
    result
}
loadone <- function(name="unif200", idx=1, expr=NULL, ...) loadall(name=name,expr=expr,idx=idx,...)[[1]]
RMSE <- function(name, ...) loadall(name,names=c("kevala","ts", "ts0"),...)

results <- data.frame(gamma1=numeric(), gamma2=numeric(), gamma3=numeric(), dist=character(),
                  rmse.kevala=numeric(),
                  rmse.ts=numeric(),
                  rmse.ts0=numeric(),
                  sd.rmse.kevala=numeric(),
                  sd.rmse.ts=numeric(),
                  sd.rmse.ts0=numeric(),
                  p=numeric())
coln <- names(results)
for (name in c("norm200","unifunif200","mix200")) {
    rmse <- as.data.frame(RMSE(name))
    row <- data.frame(1,0,1,name,
                      mean(rmse$kevala),
                      mean(rmse$ts),
                      mean(rmse$ts0),
                      sd(rmse$kevala),                      
                      sd(rmse$ts),                      
                      sd(rmse$ts0),
                      mean(rmse$kevala>=rmse$ts0))
    names(row) <- coln
    results <- rbind(results, row)
}
for (name in c("normb200","unifunifb200","mixb200")) {
    rmse <- as.data.frame(RMSE(name))
    row <- data.frame(0,1,3,name,
                      mean(rmse$kevala),
                      mean(rmse$ts),
                      mean(rmse$ts0),
                      sd(rmse$kevala),                      
                      sd(rmse$ts),                      
                      sd(rmse$ts0),
                      mean(rmse$kevala>=rmse$ts0))
    names(row) <- coln
    results <- rbind(results, row)
}

myr <- function(x,digits=3) format(round(x*10^digits,0)/10^digits,digits=digits)
res <- results[,c(1:4,c(5,7))]
for (i in 5:6) res[,i] <- myr(res[,i])
colnames(res) <- c("\\(\\bm{\\gamma_1}\\)",
                   "\\(\\bm{\\gamma_2}\\)",
                   "\\(\\bm{\\gamma_3}\\)",
                   "Model","RMSE(Kevala)","RMSE(TS)")
res$Model <- rep(c("\\(\\zeta_1\\sim\\mathcal{N}(0,1), \\zeta_1\\sim\\mathcal{N}(0,1)\\)",
                   "\\(\\zeta_1\\sim\\mathcal{U}(-6,6), \\zeta_1\\sim\\mathcal{U}(-6,6)\\)",
                   "\\(\\zeta_1\\sim\\mathcal{GM}(-4,4,0.5), \\zeta_2\\sim\\mathcal{N}(0,1)\\)"),
                 2)
a1 <- kable(res, row.names=FALSE, booktabs = TRUE,
            linesep=c('', '', '\\addlinespace'),
            caption = "My table", escape = FALSE) %>% 
    kable_styling(latex_options = "hold_position") %>%
    row_spec(0, bold=TRUE) %>% 
    ##column_spec(0, bold=TRUE) %>%
    ##collapse_rows(columns = 1)
    identity
cat(a1,"\n")


##################################################

name <- "normb200"
rmse <- RMSE(name)
MASS::eqscplot(rmse, ylab="TS RMSE", xlab="Kevala RMSE"); abline(a=0,b=1)
y1 <- rmse[,1]
y2 <- rmse[,2]
p1 <- hist(y1)
p2 <- hist(y2)
plot(p1, col=Col("darkblue", 0.4), xlim=range(c(y1,y2)))
plot(p2, col=Col("darkred",  0.4), xlim=range(c(y1,y2)), add=TRUE)
summary(rmse)
mean(rmse[,"ts0"]<=rmse[,"kevala"])

res <- loadone(name, idx=1)

pdf(file="res.pdf")
par(mfrow=c(2,1));
with(res, plot(true ~ u1, data=res$pred, type="n"))
dummy <- loadall(name,
                 lines(ts0 ~ u1, data=res$pred,
                       col=Col("darkblue", 0.2)))
Mean <- rowMeans(t(loadall(name, pred$ts0)))
with(res, lines(Mean ~ u1, data=res$pred, col="darkblue", lty=2,lwd=3))
with(res, lines(true ~ u1, data=res$pred, col="red", lty=3,lwd=3))

with(res, plot(true ~ u1, data=res$pred, type="n"))
dummy <- loadall(name,
                 lines(kevala ~ u1, data=res$pred,
                       col=Col("darkblue", 0.2)))
Mean <- rowMeans(t(loadall(name, pred$kevala)))
with(res, lines(Mean ~ u1, data=res$pred, col="darkblue", lty=2,lwd=3))
with(res, lines(true ~ u1, data=res$pred, col="red", lty=3,lwd=3))
dev.off()


