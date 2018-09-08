momest <- function(object,latent.var,p=coef(object),...) {
    if (missing(latent.var)) stop("specify name of latent variable of 1st measurement error model")
    y <- children(object, latent.var)
    lr <- lava:::lisrel(object, p)
    Lambda <- lr$Lambda[y,latent.var]
    nu <- p[y]
    Y <- model.frame(object)[y]
    x2 <- Y[,2] - (Lambda[2]*Y[,1]+nu[2])
    x3 <- Y[,3] - (Lambda[3]*Y[,1]+nu[3])
    rvar <- diag(lr$Theta)[y]
    ## response-variable
    y <- cbind(x2^3, x3^3, x2^2*x3, x3^2*x2,
               x2^4-6*Lambda[2]^2*rvar[1]*rvar[2],
               x3^4-6*Lambda[3]^2*rvar[1]*rvar[3],
               x2^3*x3-3*Lambda[2]*Lambda[3]*rvar[1]*rvar[2],
               x3^3*x2-3*Lambda[2]*Lambda[3]*rvar[1]*rvar[3])
    y <- colMeans(y)
    ## Design matrix
    B <- matrix(0,nrow=8,ncol=6)
    B[1,1:2] <- c(-Lambda[2]^3,1)
    B[2,c(1,3)] <- c(-Lambda[3]^3,1)
    B[3,1] <- -Lambda[3]*Lambda[2]^2
    B[4,1] <- -Lambda[2]*Lambda[3]^2
    B[5,4:5] <- c(Lambda[2]^4, 1)
    B[6,c(4,6)] <- c(Lambda[3]^4, 1)
    B[7,4] <- Lambda[2]^3*Lambda[3]
    B[8,4] <- Lambda[3]^3*Lambda[2]
    me <- solve(crossprod(B))%*%t(B)%*%y    
    ## Målefejl, epsilon, på measurements Y
    nom.coef <- Lambda^2/rvar
    denom <- sum(Lambda^2/rvar)
    aa <- nom.coef/denom
    ## r1 = a1*epsilon1 + a2*epsilon2 + a3*epsilon3
    ## E(r1^3) = a1^3*E(e1^3) + a2^3*E(e2^3) + a3^3*E(e3^3)
    e3 <- sum(aa^3*me[1:3])
    e4 <- sum(aa^4*me[4:6]) +
        6*aa[1]^2*aa[2]^2*rvar[1]*rvar[2] +
        6*aa[1]^2*aa[3]^2*rvar[1]*rvar[3] +
        6*aa[2]^2*aa[3]^2*rvar[2]*rvar[3]
    c(e3,e4)
}


BartlettScore <- function(object, y, newdata, p=coef(object)) { ## MLE, maximize eta in P(Y|eta; theta)
    ## y = Xb + Lambda eta + epsilon, var(e) = Psi
    lr <- lava:::lisrel(object, p)
    Lambda <- lr$Lambda
    Psi <- lr$Theta
    if (missing(y)) y <- model.frame(object)[,endogenous(object)]
    mu <- rbind(lr$mu[endogenous(object)])%x%cbind(rep(1,nrow(y)))
    if (ncol(lr$K)>0) {
        if (!missing(newdata)) {
            X <- newdata[,exogenous(object),drop=FALSE]        
        } else {
            X <- model.frame(object)[,exogenous(object),drop=FALSE]
        }
        mu <- mu + as.matrix(X)%*%(t(lr$K) + t(lr$Gamma)%*%t(Lambda))
    }
    P <- Inverse(t(Lambda)%*%Inverse(Psi)%*%Lambda)%*%t(Lambda)%*%Inverse(Psi)
    as.matrix(y-mu)%*%t(P)    
}

Wall <- function(m1,m2,data,std.err=FALSE,robust=FALSE) {
    f1 <- latent(m1)
    f2 <- latent(m2)
    m0 <- covariance(merge(m1,m2), f1,f2)
    e0 <- estimate(m0,data)
    p <- coef(e0)
    ff <- function(p) {
        lr <- lava:::lisrel(e0, p)
        Lambda <- lr$Lambda
        Psi <- lr$Theta
        rvar <- Inverse(t(Lambda)%*%Inverse(Psi)%*%Lambda)
        uhat <- BartlettScore(e0,p=p)
        u1 <- uhat[,f1]
        u2 <- uhat[,f2]
        v1 <- rvar[f1,f1]
        if (!robust) {
            e2 <- v1
            e3  <- 0     # E(epsilon^3)
            e4 <- 3*v1^2 # E(epsilon^4)
        } else {
            e2 <- v1
            mm <- momest(e0,latent.var=f1)
            e3 <- mm[1]
            e4 <- mm[2]
        }
        Mi <- cbind(1, u1, u1^2-v1,
                    u1, u1^2-v1, u1^3-3*u1*v1 - e3,
                    u1^2-v1, u1^3-3*u1*v1 - e3,
                    u1^4 -e4 -6*v1*u1^2 +6*v1^2 - 4*u1*e3)
                    ##u1^4-6*u1^2*v1 + 6*3*v1^3-3*v1^2)        
        M <- matrix(colMeans(Mi),ncol=3,byrow=TRUE)
        # m = u2 & u1*u2 & u1^2*u2-v1*u2
        m <- colMeans(cbind(u2, u1*u2, u1^2*u2-v1*u2))
        alpha <- as.vector(solve(M)%*%m)
        return(alpha)
    }
    alpha <- ff(p)
    if (!std.err) return(alpha)
    ##l <-  m - M%*%alpha
    ##Psi <- vcov(e0)
    ##summary(estimate(e0, ff))
}

