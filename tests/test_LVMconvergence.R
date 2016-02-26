rm(list=ls())

require(lava)
require(data.table)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

#### data ####
## parameters
n <- 300
rho <- rep(1,4) # 1:4 #

## model
m1 <- lvm()
regression(m1, Y1 ~ eta + epsilon1) <- list(rho[1],1)
regression(m1, Y2 ~ eta + epsilon2) <- list(rho[2],1)
regression(m1, Y3 ~ eta + epsilon3) <- list(rho[3],1)
regression(m1, Y4 ~ eta + epsilon4) <- list(rho[4],1)
intercept(m1 , ~ Y1+Y2) <- 1
intercept(m1 , ~ Y3+Y4) <- 1.5

mY.fit <- lvm()
regression(mY.fit) <- c(Y1,Y2,Y3,Y4) ~ eta
latent(mY.fit) <- ~eta

mX.fit <- lvm()
regression(mX.fit) <- c(X1,X2,X3,X4) ~ eta
latent(mX.fit) <- ~eta

## simulation 
set.seed(10)
sdX <- 1:4# rep(1,4)#

ls.data <- list()
for(iter in 1:3){
  ls.data[[iter]] <- data.table(sim(m1, n), id = 1:n)
  ls.data[[iter]][, id.char := as.character(id)]
  ls.data[[iter]][, X1 := Y1 + rnorm(.N, sd = sdX[1])]
  ls.data[[iter]][, X2 := Y2 + rnorm(.N, sd = sdX[2])]
  ls.data[[iter]][, X3 := Y3 + rnorm(.N, sd = sdX[3])]
  ls.data[[iter]][, X4 := Y4 + rnorm(.N, sd = sdX[4])]
}


resY.cv <- estimate(mY.fit,data.frame(ls.data[[1]]), control = list(constrain = TRUE))

#### CV
res0.cv <- estimate(mX.fit,data.frame(ls.data[[1]]), control = list(constrain = TRUE))
resOptimx.cv <- estimate(mX.fit,data.frame(ls.data[[1]]), control = list(constrain = TRUE),
                      method = "optimx1")
range(coef(res0.cv) - coef(resOptimx.cv))


res <- cvCheck(mX.fit, data.frame(ls.data[[1]]), n.init = 100, control = list(constrain = TRUE), verbose = TRUE)
summary(res)

#### cv issue
res0.ncv <- estimate(mX.fit, data.frame(ls.data[[3]]), control = list(constrain = TRUE))
resOptimx.ncv <- estimate(mX.fit, data.frame(ls.data[[3]]), control = list(constrain = TRUE),
                      method = "optimx1")
range(coef(res0.ncv) - coef(resOptimx.ncv))

cat(res0.ncv$opt$message)
cat(resOptimx.ncv$opt$message)

resBobyqa.ncv <- estimate(mX.fit, data.frame(ls.data[[3]]), control = list(constrain = TRUE, optimx.method = "bobyqa"),
                           method = "optimx1")
cat(resBobyqa.ncv$opt$message)
coef(resBobyqa.ncv)


check0 <- cvCheck(mX.fit, data.frame(ls.data[[3]]), n.init = 5, control = list(constrain = TRUE), verbose = TRUE)
summary(check0)
#             logLik          X2        X3        X4          eta       X2~eta       X3~eta     X4~eta      X1,X1     X2,X2        X3,X3        X4,X4      eta,eta
# range     2.255192 894.3781462 28.689168 148.02886 0.0000077309  963.5547560  30.90816048  159.47844 0.06529251 6.5803667  0.006193642 1.825386e+01 6.532402e-02
# min   -3015.541231  -0.0853617  1.831714 -14.24216 0.9281991484 -962.6619713 -31.00804747 -142.69702 4.17011650 0.1003025 12.147641607 6.491422e-79 7.157222e-06
# max   -3013.286039 894.2927845 30.520882 133.78670 0.9282068793    0.8927847  -0.09988699   16.78142 4.23540901 6.6806691 12.153835249 1.825386e+01 6.533118e-02

checkBobyqa <- cvCheck(m1.fit, data.frame(ls.data[[3]]), n.init = 5, control = list(constrain = TRUE, optimx.method = "bobyqa"), verbose = TRUE, method = "optimx1")
summary(checkBobyqa)
#           logLik        X2        X3        X4       eta     X2~eta      X3~eta    X4~eta     X1,X1     X2,X2       X3,X3     X4,X4      eta,eta
# range   137.2003  4.976931 3.1832348  353.8751  2.787073  3.7134999  2.03802343  438.5005  6.244464 0.0753927  0.07354176  9.587408 1.015136e-01
# min   -3150.5263 -3.631533 0.6394304 -104.2175 -1.546835 -2.8441585 -0.01307644 -225.7907  4.136221 6.6581975 12.08433982  6.255534 5.018685e-05
# max   -3013.3260  1.345398 3.8226652  249.6576  1.240238  0.8693413  2.02494698  212.7098 10.380685 6.7335902 12.15788158 15.842942 1.015638e-01

## centered version
data_scaled <- scale(as.data.frame(ls.data[[3]])[,paste0("X",1:4)])
res0.cncv <- estimate(mX.fit, data_scaled, control = list(constrain = TRUE))
resOptimx.cncv <- estimate(mX.fit, data_scaled, control = list(constrain = TRUE),
                      method = "optimx1")
range(coef(res0.cncv) - coef(resOptimx.cncv))

cat(res0.cncv$opt$message)
cat(resOptimx.cncv$opt$message)

resBobyqa.cncv <- estimate(mX.fit, data_scaled, control = list(constrain = TRUE, optimx.method = "bobyqa"),
                           method = "optimx1")
cat(resBobyqa.cncv$opt$message)
coef(resBobyqa.cncv)



