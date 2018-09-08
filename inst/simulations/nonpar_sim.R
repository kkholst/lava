library(lava)
matlabcmd <- "/usr/local/MATLAB/current/bin/matlab -nodisplay -r"
system("mkdir -p save")

m1 <- lvm( x1+x2+x3 ~ u1, latent= ~u1)
m2 <- lvm( y1+y2+y3 ~ u2, latent= ~u2)
m0 <- functional(merge(m1,m2), u2~u1, f=function(x) sin(x)+x)

onerun <- function(seed,...) {
    d <- sim(m,n=n,seed=seed)
    write.table(t(subset(d, selec=manifest(m))), file="save/d.csv", col.names=FALSE, row.names=FALSE, sep=",")
    write.table(t(subset(d, selec=latent(m))), file="save/du.csv", col.names=FALSE, row.names=FALSE, sep=",")
    ## Two-stage estimator, Holst-Budtz-Jørgensen 2018:
    system.time(val <- twostageCV(m1,m2,data=d,std.err=FALSE,mc.cores=2)) ## 5-fold CV
    system.time(val0 <- twostageCV(m1,m2,data=d,std.err=FALSE,mc.cores=2,nfolds=0)) ## Split in two (one for training), same tuning procedure used in Kevala 2017
    ## Kevala 2017:
    system(paste0(matlabcmd, ' "try; run(\'simkevala.m\'); catch; end; quit;"'), ignore.stdout=TRUE)
    nname <- paste0("save/kevala_",name,"_",seed,".csv")
    system(paste0("mv save/dkevala.csv ",nname))
    fit <- t(read.csv(nname,header=FALSE))
    dimnames(fit) <- list(NULL, c("u1","kevala","u1hat","u2hat"))
    
    pred <- predict(val$model2,newdata=data.frame(u1=fit[,1])); colnames(pred) <- "ts"
    pred <- cbind(fit[,1:2],pred)
    ord <- order(pred[,1])
    pred <- pred[ord,]
    pred0 <- predict(val0$model2,newdata=data.frame(u1=fit[,1])); colnames(pred0) <- "ts0"
    pred0 <- pred0[ord,"ts0"]    
    pred <- cbind(pred,ts0=pred0)
    pred <- as.data.frame(cbind(pred, true=functional(m,u2~u1)[[1]](pred[,1])))
    val$cv$fit <- NULL
    val0$cv$fit <- NULL
    res <- list(
        data=subset(d,select=c(u1,u2)),
        pred=pred, cv=list(val$cv, val0$cv),
        kevala.rmse=with(pred, mean((kevala-true)^2)),
        ts.rmse=with(pred, mean((ts-true)^2)),
        ts0.rmse=with(pred, mean((ts0-true)^2)))
    if (interactive()) {
        plot(u2 ~ u1, data=res$data)
        lines(kevala ~ u1, data=res$pred, col="red", lwd=2)
        lines(ts ~ u1, data=res$pred, col="blue", lwd=2)
        lines(true ~ u1, data=res$pred, col="gray", lty=2, lwd=3)
    }
    save(res, file=paste0("save/res_",name,"_",seed,".rda"))
    res
}
R <- seq(500)

n <- 200
m <- m0
distribution(m, ~u1) <- uniform.lvm(-6,6)
name <- paste0("unif",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun0(i),silent=TRUE)
}

m <- m0
n <- 200
name <- paste0("norm",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i),silent=TRUE)
}

n <- 200
m <- m0
distribution(m, ~u1) <- uniform.lvm(-6,6)
distribution(m, ~u2) <- uniform.lvm(-6,6)
name <- paste0("unifunif",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i),silent=TRUE)
}

n <- 200
m <- m0
distribution(m, ~u1) <- GM2.lvm()
name <- paste0("mix",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i), silent=TRUE)
}

n <- 200
m <- m0
distribution(m, ~u1) <- GM2.lvm()
distribution(m, ~u2) <- GM2.lvm()
name <- paste0("mixmix",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i), silent=TRUE)
}


##################################################


m1 <- lvm( x1+x2+x3 ~ u1, latent= ~u1)
m2 <- lvm( y1+y2+y3 ~ u2, latent= ~u2)
m0 <- functional(merge(m1,m2), u2~u1, f=function(x) sin(3*x)+x^2)

m <- m0
n <- 200
name <- paste0("normb",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i),silent=TRUE)
}

n <- 200
m <- m0
distribution(m, ~u1) <- uniform.lvm(-6,6)
distribution(m, ~u2) <- uniform.lvm(-6,6)
name <- paste0("unifunifb",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i),silent=TRUE)
}

n <- 200
m <- m0
distribution(m, ~u1) <- GM2.lvm()
name <- paste0("mixb",n)
for (i in R) {
    cat("################################### ", i,"\n")
    try(val <- onerun(i), silent=TRUE)
}
