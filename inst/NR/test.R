library(lava)
source("inst/NR/FCT.R") # derivatives of the log likelihood
# library(butils)
# butils:::package.source("lava", RorderDescription = FALSE)

## simulation
set.seed(10)

mSim <- lvm(Y~X1)
d <- sim(mSim,1e2)
eSim <- estimate(mSim, d)


## NR
set.seed(10)

mSim <- lvm(Y~X1+X2+X3)
d <- sim(mSim,1e2)
eSim <- estimate(mSim, d, control = list(iter.max = 2))
m.1 <- mSim

eT <- estimate(m.1, d, control = list(trace = 1))

e0 <- estimate(m.1, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 3))
e1 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 3, gamma = 1))
# iteration 1 should be impossible since the variance parameter is negative

#### pb with the computation of the log-likelihood: accept negative variance parameter ####
logLik(m.1, data = d, p = c(0,0,0,0,1))
logLik(m.1, data = d, p = c(0,0,0,0,-1))
logLik(m.1, data = d, p = c(0,0,0,0,-10))
score(m.1, data = d, p = c(0,0,0,0,-10))

e0 <- estimate(m.1, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 30, constrain = TRUE))
e2 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 30, gamma = 1, constrain = TRUE))
# iteration 1 is not optimal since it leads to an increase in objective function

#### pb with the size of the step ####
start <- c(-0.1493, 0, 0, 0, 1.197)

regObj(c(0,0,0,0,1), Y = d$Y, X = cbind(1,d[,-1]))
G <- regFD(start, Y = d$Y, X = cbind(1,d[,-1]))
I <- regSD(start, Y = d$Y, X = cbind(1,d[,-1]))

seq_gamma <- seq(-0.05,1,0.01)
profileLV <- sapply(seq_gamma, function(g){
  iterP <- as.double(start - g* solve(I) %*% G)
  regObj(iterP, Y = d$Y, X = cbind(1,d[,-1]))
})
plot(seq_gamma,profileLV)
abline(v = 0, col = "red")
abline(h = regObj(start, Y = d$Y, X = cbind(1,d[,-1])), col = "red")
# clearly gamma = 1 is not valid because it is a much too larger step



#### change in the backtrace ####
system.time(
  e0 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 30, constrain = TRUE))
)
system.time(
  e0 <- estimate(m.1, d, estimator = "gaussian", control = list(method = "nlminb2", trace = 1, iter.max = 30, constrain = TRUE))
)
system.time(
  e2 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 30, constrain = FALSE, backtrace = "Wolfe"))
)
system.time(
  e2 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 30, constrain = TRUE, backtrace = "Wolfe"))
)
coef(e2)-coef(eT)



#### test with a LVM ####

m <- lvm(list(Y1~eta+X1, Y2~eta, Y3~eta+X1+X2, Y4~eta+Z))
regression(m) <- eta ~ X1 + X2 + X3
covariance(m) <- Y1~Y2 

plot(m)
set.seed(10)
d <- sim(m, 1e2)

eT <- estimate(m, d)

e0 <- estimate(m, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 10))
e0bis <- estimate(m, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 10, backtrace = "Wolfe"))
range(coef(e0bis)-coef(eT))

e1 <- estimate(m, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 3, constrain = TRUE))
e1bis <- estimate(m, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 3, constrain = TRUE, backtrace = "Wolfe")) # not good, negative variance parameters ???
range(coef(e1bis)-coef(eT))


logLik(mSim,data=d,p=c(0,0,1))
lava.options(itol = -Inf)
eSim <- estimate(mSim, d)
logLik(mSim,data=d,p=c(0,0,-1))

install.packages("microbenchmark")
library(microbenchmark)
A <- rWishart(n=1, Sigma=diag(10),df=10)[,,1]

Inverse(A)

microbenchmark(svd(A),eigen(A))



