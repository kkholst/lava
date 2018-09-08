require(mets)
require(lava)
require(lava.nlin)
require(lava.mixture)
library(SQUAREM)
library(doParallel)
library(parallel)

system("mkdir -p data")

##################################################

source("wall.R")

simQuad <- function(fpostfixnum=1, # seed (appended to posfix)
             fpostfixalph="", # File postfix
             n=2e2, 			# Sample size
             param=c(vy1=1,vy2=1,vy3=1,	# Variance of residual terms of measurements (model 1)
                     vz1=1,vz2=1,vz3=1,	# Variance of residual terms of measurements (model 2)
                     bx1=-1,bx2=1,	# Covariate e
                     gam1=1,gam2=0.5,	# f2 = gam1*f1 + gam2*f1^2 + zeta2 (f1, f2 latent variables)
                     m2=1,v2=1,		# Mean and variance of residual term zeta2
                     m1a=-1,m1b=3,	# Mean of GMM for zeta1 (f1 = bx1*x + zeta1) 
                     v1a=1,v1b=1,	# Variance of GMM for zeta1
                     p1a=1),		# Mixture parameter of the GMM
             trace=0,
             mle=TRUE,
             covariate=TRUE, # Use covariates in models
             SAS=FALSE,      # Run SAS scripts
             SAS.remote=TRUE # Run SAS remotely or locally
             ) {

    suppressPackageStartupMessages(require(mets))
    suppressPackageStartupMessages(require(lava))
    suppressPackageStartupMessages(require(lava.nlin))
    
    mysummary <- function(x,sandwich=FALSE,...) {
        logL <- tryCatch(logLik(x),error=function(x) NULL)
        if (sandwich) {
            cc <- tryCatch(estimate(x)$coefmat,error=function(x) NULL)
        } else {
            cc <- coef(x,3,...)
        }
        score <- tryCatch(score(x),error=function(x) NULL)
        list(coef=cc, logLik=logL, score=score)
    }

    ## Random seed
    set.seed(fpostfixnum)
    fpostfix <- paste0(fpostfixnum,fpostfixalph)
    
    ## Simulation model
    m <- lvm()
    regression(m, c(y1,y2,y3)~u2) <- 1
    regression(m, c(z1,z2,z3)~u1) <- 1
    latent(m) <- ~u1+u2
    variance(m,~y1+y2+y3) <- c("vy1","vy2","vy3")
    variance(m,~z1+z2+z3) <- c("vz1","vz2","vz3")
    regression(m, c(u2,u1) ~ x) <- c("bx2","bx1")
    intercept(m,endogenous(m)) <- 0
    intercept(m,~u2) <- "m2"
    variance(m,~u2) <- "v2"
    transform(m,u1sq~u1) <- function(x) x^2
    regression(m,u2~u1+u1sq) <- c("gam1","gam2")
    distribution(m,parname=c("p1a","m1a","m1b","v1a","v1b"),
                 init=c(0.5,-1,1,1,1),~u1) <-
        function(n,p1a,m1a,m1b,v1a,v1b,...) {
            y1 <- rnorm(n,m1a,v1a^0.5)
            if (p1a>=1) return(y1)
            z <- rbinom(n,1,p1a)
            y2 <- rnorm(n,m1b,v1b^0.5)
            return(z*y1+(1-z)*y2)
        }

    d <- subset(d. <- sim(m,n,p=param),
                select=c(y1,y2,y3,z1,z2,z3,x))
    
    ## Save data and parameters to csv
    dfile <- paste0("data/data_",fpostfix,".csv")
    pfile <- paste0("data/param_",fpostfix,".csv")
    write.csv(d,file=dfile,row.names=FALSE)
    write.csv(rbind(param),file=pfile,row.names=FALSE)

    ##################################################
    
    ## Linear model ignoring non-linear effect of u1 on u2
    m0 <- lvm(c(z1,z2,z3)~u1,c(y1,y2,y3)~u2,u2~u1+x,u1~x)
    if (!covariate) rmvar(m0) <- ~x
    e0 <- tryCatch(estimate(m0,data=d),
                   error=function(x) NULL)

    ##################################################
    
    ## MLE (Adaptive Gaussian Quadrature)
    p0 <- with(as.data.frame(rbind(param)),
               c(m1a,0,0,m2,0,0, ## Intercepts
                 1,1,1,1, ## Factor loadings
                 if(covariate) c(bx1,bx2),
                 gam1,gam2, ## Regression coefficients
                 log(v1a),log(v2),0,0,0,0,0,0 ## Residual variance (log)
                 ))    
    if (covariate) {
        predform <- ~x
    } else {
        predform <- ~1
    }

    e1 <- e1a <- e1b <- NULL
    if (mle) {
        nmodel <- list(measure1=~z1+z2+z3,measure2=~y1+y2+y3,pred1=predform,pred2=predform,model="nsem2")
        ## Laplace approximation, nq=1
        e1 <- tryCatch(nsem(nmodel,data=d,control=list(trace=trace,start=p0),laplace.control=list(nq=0,lambda=0.3,eb0=FALSE)),
                       error=function(x) NULL)
        p1 <- p0
        if (!is.null(e1)) p1 <- coef(e1)
    ## Non-linear model, MLE AQG with nq=3, reuse empirical bayes estimates
        e1a <- tryCatch(nsem(nmodel,data=d,control=list(trace=trace,start=p0),laplace.control=list(nq=3,lambda=0.3,eb0=TRUE)),
                    error=function(x) NULL)
        ## Non-linear model, MLE AQG with nq=9
        e1d <- tryCatch(nsem(nmodel,data=d,control=list(trace=trace,start=p1),laplace.control=list(nq=9,lambda=0.3,eb0=FALSE)),
                        error=function(x) NULL)
    }
    
    ##################################################

    ## Two stage model
    m1 <- lvm(c(z1,z2,z3)~u1,u1~x,latent=~u1)
    if (!covariate) {
        rmvar(m1) <- ~x
    }
    M1 <- estimate(m1,d)
    m2 <- lvm(c(y1,y2,y3)~u2,u2~x,latent=~u2)
    if (!covariate) {
        rmvar(m2) <- ~x
    }    
    nonlinear(m2,type="quadratic") <- u2~u1
    ts1 <- tryCatch(lava::twostage(M1,m2,data=d),
                    error=function(x) NULL)

    ## Two stage mixture
    m1a <- baptize(fixsome(m1))
    intercept(m1a,latent(m1a)) <- NA
    p2 <- with(as.data.frame(rbind(param)),
               c(0,0,m1a,m1b, ## Intercepts
                 1,1,         ## Factor loadings
                 if (covariate) bx1,         ## Regression coefficients
                 1,1,1,v1a,   ## Residual variance (log)
                 abs(p1a-0.1)     ## Class 1 probability
                 ))
    ## mm <- mixture(list(m1a,m1a),data=d,control=list(K=1,constrain=TRUE,start=p2,trace=1))    
    M2  <- tryCatch(mixture(list(m1a,m1a),data=d,control=list(start=p2,trace=trace)),
                    error=function(x) NULL)
    ts2 <- tryCatch(lava::twostage(M2,m2,data=d),
                    error=function(x) NULL)

    ##################################################
    
    ## IV estimators
    d2 <- transform(d,z12=z1^2,z22=z2^2,z32=z3^2)
    ## Single measure
    miv1 <- lvm(c(z1,z2,z3)~f1, c(z12,z22,z32)~f2, f1~x, f2~x,
                y1~f1+f2+x)
    if (!covariate) {
        rmvar(miv1) <- ~x
    }
    eiv1 <- tryCatch(estimate(miv1,d2,estimator="iv0"),error=function(x) NULL)
    eiv1.rob <- tryCatch(estimate(miv1,d2,estimator="iv"),error=function(x) NULL)

    ## Multiple measures (same?)
    miv2 <- lvm(c(z1,z2,z3)~f1, c(z12,z22,z32)~f2, f1~x, f2~x,
           c(y1,y2,y3) ~eta,
           eta~f1+f2+x)
    if (!covariate) {
        rmvar(miv2) <- ~x
    }

    eiv2 <- tryCatch(estimate(miv2,d2,estimator="iv"),error=function(x) NULL)
    ## Pseudo-likelihood
    eiv2.pmle <- tryCatch(estimate(miv2,d2),error=function(x) NULL)
    eiv2.pmle.sandwich <- tryCatch(estimate(eiv2.pmle),error=function(x) NULL)

##################################################
    
    wall <- tryCatch(Wall(m1,m2,data=d), error=function(x) NULL)
    wall.rob <- tryCatch(Wall(m1,m2,data=d,robust=TRUE), error=function(x) NULL)
    if (length(wall)==3L) {
        wall <- cbind(Estimate=wall, "Std. Error"=NA)
        wall.rob <- cbind(Estimate=wall.rob, "Std. Error"=NA)
        rownames(wall) <- rownames(wall.rob) <- c("(Intercept)","f1","f1^2")
    }
    
##################################################

    sasIV <- sasNSEM <- sasNSEMMIX <- NULL
    if (SAS) {
        ## SAS (nonlinear sem, nonlinear mixture sem, IV)
        sascmd <- paste0("make ","sas",if (SAS.remote) "r"," postfix=",fpostfix)
        try(system(sascmd,intern=TRUE))
        
        fsasname <- paste0("data/sas_",fpostfix,"_")
        sasIV <- tryCatch(lava:::getSAS(paste0(fsasname,"iv.csv")),error=function(x) NULL)
        sasNSEM <- tryCatch(lava:::getSAS(paste0(fsasname,"nonlin.csv")),error=function(x) NULL)
        sasNSEMMIX <- tryCatch(lava:::getSAS(paste0(fsasname,"mixnonlin.csv")),error=function(x) NULL)
    }         

##################################################

    res <- list(linear=mysummary(e0),
                
                laplace=if (mle) mysummary(e1) else NULL,
                agq3=if (mle) mysummary(e1a) else NULL,
                agq9=if (mle) mysummary(e1d) else NULL,
                twostage=coef(ts1,3),
                twostage.naive=coef(ts1$naive,3),
                twostage.naive.robust=ts1$naive.robust$coefmat,
                twostage.mixture=coef(ts2,3),
                model1=mysummary(M1,sandwich=TRUE),
                model1.mixture=mysummary(M2,sandwich=TRUE),
                iv1=coef(eiv1,3),
                iv1.rob=coef(eiv1.rob,3),
                iv2=coef(eiv2,3),
                iv.pmle=coef(eiv2.pmle,3),
                iv.pmle.robust=eiv2.pmle.sandwich$coefmat,
                wall=wall,
                wall.rob=wall.rob,
                sas.iv=sasIV,
                sas.nsem=sasNSEM,
                sas.mixture=sasNSEMMIX)

    save(res,file=paste0("data/res_",fpostfix,".rda"))

    invisible(res)
}
