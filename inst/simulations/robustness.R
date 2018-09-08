library(mets)

##################################################
## Exponential
##################################################

beta <- 0.3
f <- function(x) beta*exp(x)
m <- lvm(x1+x2+x3~eta1, y1+y2+y3~eta2, latent=~eta1+eta2)
variance(m,~eta2+eta1) <- c(0.2,0.5)
functional(m, eta2~eta1) <- f
d <- sim(m,500,seed=1,latent=TRUE)
plot(eta2~eta1,d)

m1 <- lvm(x1+x2+x3~eta1,latent=~eta1)
m2 <- lvm(y1+y2+y3~eta2,latent=~eta2)
nonlinear(m2,type="exp") <- eta2~eta1
a <- twostage(m1,m2, data=d)

d2 <- transform(d,x12=exp(x1),x22=exp(x2),x32=exp(x3))
miv <- lvm(c(x1,x2,x3)~f1, c(x12,x22,x32)~f2, f1~1, f2~1,
          y1~f1+f2)
eiv <- tryCatch(estimate(miv,d2,estimator="iv"),error=function(x) NULL)

onerun <- function(...,n=500) {
    d <- sim(m,500,latent=FALSE)
    a <- twostage(m1,m2, data=d)   
    res <- estimate(a,robust=FALSE,keep="eta2~eta1",regex=TRUE)    
    d2 <- transform(d,x12=exp(x1),x22=exp(x2),x32=exp(x3))
    eiv <- tryCatch(estimate(miv,d2,estimator="iv"),error=function(x) NULL)
    res2 <- estimate(eiv,robust=FALSE,keep=c("y1~f1","y1~f2"))
    c(res$coefmat[,1],res2$coefmat[,1],
      res$coefmat[,2],res2$coefmat[,2])
}

set.seed(1)
val <- sim(onerun,1000,mc.cores=min(detectCores()-1,10),messages=1,args=list(n=500))
save(val, file="data/exp1.rda")


##################################################
## Quadratic, with uniform distributions
##################################################

f <- function(x) x+0.5*x*x
m <- lvm(x1+x2+x3~eta1, y1+y2+y3~eta2, latent=~eta1+eta2)
## variance(m,~eta2+eta1) <- c(0.2,0.5)
functional(m, eta2~eta1) <- f
d <- sim(m,500,seed=1,latent=TRUE)
m1 <- lvm(x1+x2+x3~eta1,latent=~eta1)
m1 <- baptize(m1)
intercept(m1, ~x1+eta1) <- c(0,NA)
regression(m1, x1~eta1) <- 1
m2 <- lvm(y1+y2+y3~eta2,latent=~eta2)
nonlinear(m2,type="quadratic") <- eta2~eta1
a <- twostage(m1,m2, data=d)



#g0 <- estimate(m1,d)
p0 <- c(x2=0,x3=0,p1=-1,p2=1,
       "x2~eta1"=1,"x3~eta1"=1,
       "x1~~x1"=1,
       "x2~~x2"=1,
       "x3~~x3"=1,
       "eta1~~eta1"=0.5,
       "pr1"=0.5)
## g <- mixture(list(m1,m1),data=d,control=list(trace=0,start=p0))
## twostage(g,m2,data=d)
K <- 3
p3 <- c(x2=0,x3=0,
       seq(-1.7,1.7,length.out=K),
       "x2~eta1"=1,
       "x3~eta1"=1,
       "x1~~x1"=1,
       "x2~~x2"=1,
       "x3~~x3"=1,
       "eta1~~eta1"=0.5,
       rep(1/K, K-1))
g3 <- mixture(rep(list(m1),3),data=d,control=list(trace=0,start=p3))
K <- 4
p4 <- c(x2=0,x3=0,
       seq(-1.7,1.7,length.out=K),
       "x2~eta1"=1,
       "x3~eta1"=1,
       "x1~~x1"=1,
       "x2~~x2"=1,
       "x3~~x3"=1,
       "eta1~~eta1"=0.5,
       rep(1/K, K-1))


m0 <- lvm(x~1)
covariance(m0,~x) <- "v"
distribution(m0,~x) <- uniform.lvm()
d0 <- sim(m0,500)
E1 <- estimate(m0,d0)
K <- 2; E2 <- mixture(rep(list(m0),K),data=d0,control=list(trace=1,start=c(seq(-1.7,1.7,length.out=K),v=0.5,rep(1/K,K-1))))
K <- 3; E3 <- mixture(rep(list(m0),K),data=d0,control=list(trace=1,start=c(seq(-1.7,1.7,length.out=K),v=0.5,rep(1/K,K-1))))
K <- 4; E4 <- mixture(rep(list(m0),K),data=d0,control=list(trace=1,start=c(seq(-1.7,1.7,length.out=K),v=0.5,rep(1/K,K-1))))
K <- 5; E5 <- mixture(rep(list(m0),K),data=d0,control=list(trace=1,start=c(seq(-1.7,1.7,length.out=K),v=0.5,rep(1/K,K-1))))
AIC(E1,E2,E3,E4,E5)
density(d0$x)

onerun <- function(...,n=500) {
    d <- sim(m,500,latent=FALSE)
    a <- twostage(m1,m2, data=d)
    res <- estimate(a,robust=FALSE,keep="eta2~eta1",regex=TRUE)
    g2 <- mixture(list(m1,m1),data=d,control=list(trace=0,start=p0))
    b2 <- twostage(g2,m2, data=d)
    res2 <- estimate(b2,robust=FALSE,keep="eta2~eta1",regex=TRUE)
    g3 <- mixture(rep(list(m1),3),data=d,control=list(trace=0,start=p3))
    b3 <- twostage(g3,m2,data=d)
    res3 <- estimate(b3,robust=FALSE,keep="eta2~eta1",regex=TRUE)
#    g4 <- mixture(rep(list(m1),4),data=d,control=list(trace=0,start=p4))
#    b4 <- twostage(g4,m2,data=d)
#    res4 <- estimate(b4,robust=FALSE,keep="eta2~eta1",regex=TRUE)    
    c(res$coefmat[,1],res2$coefmat[,1],res3$coefmat[,1],#res4$coefmat[,1],
      res$coefmat[,2],res2$coefmat[,2],res3$coefmat[,2])#res4$coefmat[,2])
}

m1 <- baptize(m1)
intercept(m1, ~x1+eta1) <- c(0,NA)
regression(m1, x1~eta1) <- 1
distribution(m,~x1+x2+x3) <- uniform.lvm()
set.seed(1)

#val <- sim(onerun,20,mc.cores=4,messages=1,args=list(n=500))
#val <- sim(val,100,mc.cores=4,messages=1,args=list(n=500))
val1 <- sim(onerun,1000,mc.cores=min(detectCores()-1,20),messages=1,args=list(n=500))
save(val1, file="data/uniform1.rda")


##################################################

m1 <- lvm(x1+x2+x3~eta1,latent=~eta1)
m1 <- baptize(m1)
intercept(m1, ~x1+eta1) <- c(0,NA)
regression(m1, x1~eta1) <- 1
distribution(m,~eta1) <- uniform.lvm()
set.seed(1)

val2 <- sim(onerun,1000,mc.cores=min(detectCores()-1,20),messages=1,args=list(n=500))
save(val2, file="data/uniform2.rda")


