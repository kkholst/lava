source("nlinsem.R")

seqR <- seq(1,1000)
n <- 1e3
mc.cores <- min(detectCores()-1,38)

onerun <- function(x,...) {
    message("*** Iter ",x)
    val <- tryCatch(
        simQuad(x,fpostfixalph,
                param=param,n=n,
                ...),
        error=function(x) NULL)
}


fpostfixalph <- "_normalc"
param <- c(vy1=1,vy2=1,vy3=1,
           vz1=1,vz2=1,vz3=1,
           bx1=0,bx2=0,
           gam1=1,gam2=0.5,
           m2=1,v2=1,
           m1a=-3,m1b=3,
           v1a=1,v1b=1,
           p1a=1)
res <- mclapply(seqR,onerun, mc.cores=mc.cores,covariate=FALSE)

fpostfixalph <- "_mixture2c"
param <- c(vy1=1,vy2=1,vy3=1,
           vz1=1,vz2=1,vz3=1,
           bx1=0,bx2=0,
           gam1=1,gam2=0.5,
           m2=1,v2=1,
           m1a=0,m1b=3,
           v1a=1,v1b=1,
           p1a=0.25)
res <- mclapply(seqR,onerun, mc.cores=mc.cores,covariate=FALSE,mle=FALSE)


################################################################################


fpostfixalph <- "_normalxc"
param <- c(vy1=1,vy2=1,vy3=1,
           vz1=1,vz2=1,vz3=1,
           bx1=-1,bx2=1,
           gam1=1,gam2=0.5,
           m2=1,v2=1,
           m1a=-3,m1b=3,
           v1a=1,v1b=1,
           p1a=1)
res <- mclapply(seqR,onerun, mc.cores=mc.cores,covariate=TRUE)

fpostfixalph <- "_mixture2xc"
param <- c(vy1=1,vy2=1,vy3=1,
           vz1=1,vz2=1,vz3=1,
           bx1=-1,bx2=1,
           gam1=1,gam2=0.5,
           m2=1,v2=1,
           m1a=0,m1b=3,
           v1a=1,v1b=1,
           p1a=0.25)
res <- mclapply(seqR,onerun, mc.cores=mc.cores,covariate=TRUE,mle=FALSE)



