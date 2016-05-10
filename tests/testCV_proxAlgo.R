rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)
library(lava)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

coef2.penalized <- function(x, iter_lambda){
  
  if(is.list(x)){
    
    if(!missing(iter_lambda)){
      x <- x[[iter_lambda]]
    }else{
      res <- lapply(x, function(model){
        c(model@lambda1,
          model@lambda2,
          model@unpenalized, 
          model@penalized, 
          model@nuisance$sigma2)
      })
      
      Mres <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
      colnames(Mres) <- c("lambda1","lambda2",names(x[[1]]@unpenalized),names(x[[1]]@penalized),"sigma2")
      return(Mres)
    }
    
  } 
  
  coef <- c(x@unpenalized, 
            x@penalized, 
            x@nuisance$sigma2)
  return(coef)
}

validLVM <- function(x){
  library(penalized)
  
  name.outcome <- endogenous(x)
  lambda <- x$penalty$lambda1*coef(x)[paste(name.outcome,name.outcome,sep = ",")]
  resPenalized <- penalized(as.formula(paste0(name.outcome,"~.")), lambda1 = lambda, data = x$data$model.frame)
  diffCoef <- coef(x) - coef2.penalized(resPenalized)
  
  return(diffCoef)
}

#### 0- problematic examples ####

#### data
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

#### models
lvm.model <- lvm(formula.lvm)
elvm.model <- estimate(lvm.model, df.data)
plvm.model <- penalize(lvm.model)

#### regularization path
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seqPark_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
seqParkNorm_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1/x@nuisance$sigma2}))

#### specific knot
lambda_tempo <- seqParkNorm_lambda[6]

start1 <- coef(estimate(lvm.model, data = df.data))
start2 <- coef(estimate(plvm.model, lambda1 = 1e5, data = df.data))

plvm.punctual1 <- estimate(plvm.model, data = df.data, lambda1 = lambda_tempo,
                           control = list(start = start1))

plvm.punctual2 <- estimate(plvm.model, data = df.data, lambda1 = lambda_tempo,
                           control = list(start = start2))

validLVM(plvm.punctual1)
validLVM(plvm.punctual2)
coef(plvm.punctual1) - coef(plvm.punctual2)

#### 1- Definition of the path ####

set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)

system.time(
  lvm.test0 <- estimate(plvm.model,  data = df.data)
)

system.time(
  lvm.testPath <- estimate(plvm.model,  data = df.data, 
                           regularizationPath =  TRUE,
                           control = list(constrain = FALSE, iter.max = 5000, fast = 3, eta.BT = 0.8))
)
lvm.testPath$opt$message

#### 2- at one node ####
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)
seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))
n.lambda <- length(seq_lambda)

#### only L1
df.res <- data.frame(matrix(nrow = n.lambda*10, ncol  = 6))
names(df.res) <- c("step","eta.BT","lambda","time","iteration","maxDiff")

for(iter_l in 1:length(seq_lambda)){
  cat(iter_l, " : ",penalized.PathL1[[iter_l]]@lambda1,"\n")
  
  time <- system.time(
    fit_tempo <- estimate(plvm.model,  data = df.data, 
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                             control = list(constrain = TRUE, iter.max = 1000, step = NULL))
  )
  
  df.res[(iter_l-1)*n.lambda + 1,"step"] <- "Hessian"
  df.res[(iter_l-1)*n.lambda + 1,c("lambda","time","iteration","maxDiff")] <- c(penalized.PathL1[[iter_l]]@lambda1, time[3], fit_tempo$opt$iterations, max(abs(validLVM(fit_tempo))))
  
  time <- system.time(
    fit_tempo <- estimate(plvm.model,  data = df.data, 
                          lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 1000, step = NULL, fast = 3))
  )
  
  df.res[(iter_l-1)*n.lambda + 2,"step"] <- "Hessian-Fast3"
  df.res[(iter_l-1)*n.lambda + 2,c("lambda","time","iteration","maxDiff")] <- c(penalized.PathL1[[iter_l]]@lambda1, time[3], fit_tempo$opt$iterations, max(abs(validLVM(fit_tempo))))
  
  
  time <- system.time(
    fit_tempo <- estimate(plvm.model,  data = df.data, 
                          lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 1000, step = 1, n.BT = 10, eta.BT = 0.5))
  )
  
  df.res[(iter_l-1)*n.lambda + 3,"step"] <- "BT"
  df.res[(iter_l-1)*n.lambda + 3,c("eta.BT", "lambda","time","iteration","maxDiff")] <- c(0.5, penalized.PathL1[[iter_l]]@lambda1, time[3], fit_tempo$opt$iterations, max(abs(validLVM(fit_tempo))))
  
  time <- system.time(
    fit_tempo <- estimate(plvm.model,  data = df.data, 
                          lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 1000, step = 1, n.BT = 30, eta.BT = 0.8))
  )
  
  df.res[(iter_l-1)*n.lambda + 4,"step"] <- "BT"
  df.res[(iter_l-1)*n.lambda + 4,c("eta.BT", "lambda","time","iteration","maxDiff")] <- c(0.8, penalized.PathL1[[iter_l]]@lambda1, time[3], fit_tempo$opt$iterations, max(abs(validLVM(fit_tempo))))
  
  time <- system.time(
    fit_tempo <- estimate(plvm.model,  data = df.data, 
                          lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                          control = list(constrain = TRUE, iter.max = 1000, step = 1, n.BT = 30, eta.BT = 0.8, fast = 3))
  )
  
  df.res[(iter_l-1)*n.lambda + 5,"step"] <- "BT-fast3"
  df.res[(iter_l-1)*n.lambda + 5,c("eta.BT", "lambda","time","iteration","maxDiff")] <- c(0.8, penalized.PathL1[[iter_l]]@lambda1, time[3], fit_tempo$opt$iterations, max(abs(validLVM(fit_tempo))))
  
}





# step eta.BT     lambda time iteration      maxDiff
# 1        Hessian     NA 1617.81111 0.13         2 1.065814e-14
# 2  Hessian-Fast3     NA 1617.81111 0.14         2 8.326673e-17
# 3             BT    0.5 1617.81111 0.14         2 1.492140e-13
# 4             BT    0.8 1617.81111 0.15         2 2.096101e-13
# 5       BT-fast3    0.8 1617.81111 0.14         2 3.019807e-13
# 6           <NA>     NA         NA   NA        NA           NA
# 7           <NA>     NA         NA   NA        NA           NA
# 8           <NA>     NA         NA   NA        NA           NA
# 9           <NA>     NA         NA   NA        NA           NA
# 10       Hessian     NA 1617.64933 0.13         2 1.065814e-14
# 11 Hessian-Fast3     NA 1617.64933 0.12         2 8.326673e-17
# 12            BT    0.5 1617.64933 0.13         2 1.215028e-12
# 13            BT    0.8 1617.64933 0.13         2 1.215028e-12
# 14      BT-fast3    0.8 1617.64933 0.11         2 3.019807e-13
# 15          <NA>     NA         NA   NA        NA           NA
# 16          <NA>     NA         NA   NA        NA           NA
# 17          <NA>     NA         NA   NA        NA           NA
# 18          <NA>     NA         NA   NA        NA           NA
# 19       Hessian     NA 1002.92748 2.69       292 7.098295e-08
# 20 Hessian-Fast3     NA 1002.92748 3.04       304 7.037149e-08
# 21            BT    0.5 1002.92748 0.90       112 1.509833e-06
# 22            BT    0.8 1002.92748 1.26       141 8.984192e-07
# 23      BT-fast3    0.8 1002.92748 2.49       281 2.568184e-08
# 24          <NA>     NA         NA   NA        NA           NA
# 25          <NA>     NA         NA   NA        NA           NA
# 26          <NA>     NA         NA   NA        NA           NA
# 27          <NA>     NA         NA   NA        NA           NA
# 28       Hessian     NA  606.39734 4.26       373 2.803123e-08
# 29 Hessian-Fast3     NA  606.39734 4.29       385 2.898022e-08
# 30            BT    0.5  606.39734 1.31       168 4.106298e-07
# 31            BT    0.8  606.39734 1.68       197 3.613158e-07
# 32      BT-fast3    0.8  606.39734 3.01       330 7.058585e-08
# 33          <NA>     NA         NA   NA        NA           NA
# 34          <NA>     NA         NA   NA        NA           NA
# 35          <NA>     NA         NA   NA        NA           NA
# 36          <NA>     NA         NA   NA        NA           NA
# 37       Hessian     NA  606.33670 9.72      1001 1.019603e-06
# 38 Hessian-Fast3     NA  606.33670 9.97      1001 1.073987e-06
# 39            BT    0.5  606.33670 3.85       548 8.796532e-07
# 40            BT    0.8  606.33670 4.34       585 1.217964e-06
# 41      BT-fast3    0.8  606.33670 8.69      1001 9.259773e-07
# 42          <NA>     NA         NA   NA        NA           NA
# 43          <NA>     NA         NA   NA        NA           NA
# 44          <NA>     NA         NA   NA        NA           NA
# 45          <NA>     NA         NA   NA        NA           NA
# 46       Hessian     NA   29.44781 0.42        32 5.459521e-08
# 47 Hessian-Fast3     NA   29.44781 0.52        38 7.425858e-08
# 48            BT    0.5   29.44781 0.29        17 3.771550e-08
# 49            BT    0.8   29.44781 0.31        16 3.454156e-07
# 50      BT-fast3    0.8   29.44781 0.42        20 3.933269e-09
# 51          <NA>     NA         NA   NA        NA           NA
# 52          <NA>     NA         NA   NA        NA           NA
# 53          <NA>     NA         NA   NA        NA           NA
# 54          <NA>     NA         NA   NA        NA           NA
# 55       Hessian     NA   29.44487 0.38        32 5.458809e-08
# 56 Hessian-Fast3     NA   29.44487 0.48        38 7.425106e-08
# 57            BT    0.5   29.44487 0.37        18 7.619413e-08
# 58            BT    0.8   29.44487 0.32        16 3.453575e-07
# 59      BT-fast3    0.8   29.44487 0.37        20 1.037084e-08
# 60          <NA>     NA         NA   NA        NA           NA
# 61          <NA>     NA         NA   NA        NA           NA
# 62          <NA>     NA         NA   NA        NA           NA
# 63          <NA>     NA         NA   NA        NA           NA
# 64       Hessian     NA   24.99796 0.37        32 4.967298e-08
# 65 Hessian-Fast3     NA   24.99796 0.47        39 3.771826e-08
# 66            BT    0.5   24.99796 0.23        15 9.473576e-08
# 67            BT    0.8   24.99796 0.33        18 3.377237e-08
# 68      BT-fast3    0.8   24.99796 0.35        18 1.830617e-07
# 69          <NA>     NA         NA   NA        NA           NA
# 70          <NA>     NA         NA   NA        NA           NA
# 71          <NA>     NA         NA   NA        NA           NA
# 72          <NA>     NA         NA   NA        NA           NA
# 73       Hessian     NA    0.00000 0.32        23 1.839351e-05
# 74 Hessian-Fast3     NA    0.00000 0.39        29 1.434730e-05
# 75            BT    0.5    0.00000 0.20        12 1.629022e-05
# 76            BT    0.8    0.00000 0.28        13 1.262891e-05
# 77      BT-fast3    0.8    0.00000 0.33        15 7.168167e-06

