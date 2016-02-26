#%%%%%
#%
#% Test convergence of the penalized LVM functions
#%
#%%%%%

rm(list = ls())

set.seed(10)

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)

path.lava <- "C:/Users/hpl802/Documents/GitHub/lava"
vecRfiles <- list.files(file.path(path.lava,"R"))
sapply(vecRfiles, function(x){source(file.path(path.lava,"R",x))})

#### Over the grid ####

## parametrisation
seq_n <- c(25,50,100,200,300,500)
seq_p <- c(5,10,15,20)
seq_lambda <- c(1,5,10,15)
seq_rep <- 1:5
Allbeta <- beta <- c(1:5,rep(0, max(seq_p) - 5))# rbinom(p, size = 1, prob = 1/seq(1,p/2, by = 0.5))
method <- "nlminb2" # "optimx1" #"NR"

grid <- expand.grid(n =  seq_n,
                    p =  seq_p,
                    lambda =  seq_lambda,
                    rep =  seq_rep)
n.grid <- nrow(grid)

dt.results <- data.table(matrix(as.numeric(NA), nrow = n.grid, ncol = 16))
names(dt.results) <- c("n.obs","n.param","lambda","rep",
                       "L1.zeroBeta","L1.cv","L1.iteration","L1.maxDiff","L1.minDiff",
                       "L2.cv","L2.iteration","L2.maxDiff","L2.minDiff",
                       "time.lvm","time.lvmL1","time.lvmL2")

.pb <- txtProgressBar(min = 0, max = n.grid, style = 3)


for(iter_grid in 1:n.grid){

  setTxtProgressBar(.pb, iter_grid)
  
  n <- grid[iter_grid,"n"]
  p <- grid[iter_grid,"p"]
  lambda1 <- grid[iter_grid,"lambda"]
  lambda2 <- grid[iter_grid,"lambda"] 
  beta <- Allbeta[1:p] 
  
  ## simulation
  X <- matrix(rnorm(n*p),n,p)
  Y <- rowSums( sweep(X, MARGIN = 2, STATS= beta, FUN = "*") )+ rnorm(n)
  df.data <- setNames(data.frame(Y,X),
                      c("Y", paste0("X",1:p))
  )
  
  ## formula
  formula_tempo <- as.formula(paste0("Y ~ ", paste(paste0("X",1:p),collapse="+")))
  lvm.model <- lava:::lvm( formula_tempo )

  penalty1 <- ls.fct_penalty(lvm = lvm.model,
                             lambda1 = lambda1, lambda2 = 0)
  
  penalty2 <- ls.fct_penalty(lvm = lvm.model,
                             lambda1 = 0, lambda2 = lambda2)

  ## lava
  tps1 <- system.time(
  elvm.fit <- estimate(lvm.model, df.data, 
                       control = list(constrain = TRUE))
  )
  tps2 <- system.time(
  elvm1.fit <- estimate(lvm.model, data = df.data, estimator = "penalised" ,penalisation = penalty1, # df.dataAll[,all.vars(formula_tempo)]
                        control = list(constrain = TRUE, method = method))
  )
  tps3 <- system.time(
  elvm2.fit <- estimate(lvm.model,  data = df.data, estimator = "penalised" ,penalisation = penalty2,
                        control = list(constrain = TRUE, method = method))
  )
  
  ## penalized
  penalized1.fit <- penalized:::penalized(response = df.data$Y, penalized = df.data[,all.vars(formula_tempo)[-1]], 
                                          lambda1 = lambda1 * coef(elvm1.fit)["Y,Y"], lambda2 = 0, trace = FALSE )
  
  penalized2.fit <- penalized:::penalized(response = df.data$Y, penalized = df.data[,all.vars(formula_tempo)[-1]], 
                                          lambda1 = 0, lambda2 = lambda2 * coef(elvm2.fit)["Y,Y"], trace = FALSE )
  
  
  dt.results[iter_grid, c("n.obs","n.param","lambda","rep") := .(n,p,grid[iter_grid,"lambda"],grid[iter_grid,"rep"])]
  dt.results[iter_grid, c("time.lvm","time.lvmL1","time.lvmL2") := .(tps1[3],tps2[3],tps3[3])]
  
  param_penalized1 <- c(penalized1.fit@unpenalized, penalized1.fit@penalized,  penalized1.fit@nuisance$sigma2)
  dt.results[iter_grid, L1.zeroBeta := sum( abs(penalized1.fit@penalized) < 1e-10)]
  dt.results[iter_grid, L1.cv := if("iterations" %in% names(elvm1.fit)){elvm1.fit$opt$convergence}else{NA} ]
  dt.results[iter_grid, L1.iteration := if("iterations" %in% names(elvm1.fit)){elvm1.fit$opt$iterations}else{NA} ]
  dt.results[iter_grid, L1.maxDiff := max(abs(param_penalized1-coef(elvm1.fit)))]
  dt.results[iter_grid, L1.minDiff := min(abs(param_penalized1-coef(elvm1.fit)))]
  
  param_penalized2 <- c(penalized2.fit@unpenalized, penalized2.fit@penalized,  penalized2.fit@nuisance$sigma2)
  dt.results[iter_grid, L2.cv := if("iterations" %in% names(elvm2.fit)){elvm2.fit$opt$convergence}else{NA} ]
  dt.results[iter_grid, L2.iteration := if("iterations" %in% names(elvm2.fit)){elvm2.fit$opt$iterations}else{NA} ]
  dt.results[iter_grid, L2.maxDiff := max(abs(param_penalized2-coef(elvm2.fit)))]
  dt.results[iter_grid, L2.minDiff := min(abs(param_penalized2-coef(elvm2.fit)))]
  
  
}
dt.results

save(dt.results, file = paste0("C:/Users/hpl802/Documents/Projects/LVM/CV-test/dt_results-",method,".RData"))

#### display
# dt.results[, .(L2.cv = 100 * mean(L2.cv==0)), by =c("n.obs","n.param","lambda")]
# 
# dt.results[, .(L2.cv = 100 * mean(L2.cv==0)), by =c("n.obs","n.param")]
# 
# dt.results[, .(L2.cv = 100 * mean(L2.cv==0)), by =c("n.param")]
# 
# dt.results[, .(Mismatch = 100 * mean(L2.maxDiff>1e-3)), by =c("n.obs","n.param","lambda")]
# 
# dt.results[, table(L2.cv,L2.maxDiff>1e-2)]

