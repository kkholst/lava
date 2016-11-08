library(lava)
library(reshape2)
library(plotly)
library(plot3D)
source("inst/NR/FCT.R") # derivatives of the log likelihood

## simulation
set.seed(10)

mSim <- lvm(Y~X1)
d <- sim(mSim,1e2)
eSim <- estimate(mSim, d)

## NR
set.seed(10)

mSim <- lvm(Y~X1+X2+X3)
d <- sim(mSim,1e2)
eSim <- estimate(mSim, d)
m.1 <- mSim

eT <- estimate(m.1, d, control = list(trace = 1))

e0 <- estimate(m.1, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 3))
e1 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 3, gamma = 1))
# iteration 1 should be impossible since the variance parameter is negative

#### pb with the computation of the log-likelihood: accept negative variance parameter ####
logLik(m.1, data = d, p = c(0,0,0,0,1))
logLik(m.1, data = d, p = c(0,0,0,0,-1))
logLik(m.1, data = d, p = c(0,0,0,0,-10))

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
e0 <- estimate(m.1, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 30, constrain = TRUE, backtrace = 2))
e2 <- estimate(m.1, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 30, gamma = 1, constrain = TRUE, backtrace = 2))

coef(e0)-coef(eT)



#### test with a LVM ####

m <- lvm(list(Y1~eta+X1, Y2~eta, Y3~eta+X1+X2, Y4~eta+Z))
regression(m) <- eta ~ X1 + X2 + X3
covariance(m) <- Y1~Y2 

plot(m)
set.seed(10)
d <- sim(m, 1e2)

eT <- estimate(m, d)

e0 <- estimate(m, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 10))
e0bis <- estimate(m, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 10, backtrace = 2))

e1 <- estimate(m, d, estimator = "gaussian1", control = list(method = "NR", trace = 1, iter.max = 10))
e1bis <- estimate(m, d, estimator = "gaussian", control = list(method = "NR", trace = 1, iter.max = 10, backtrace = 2))
range(coef(e1bis)-coef(eT))
