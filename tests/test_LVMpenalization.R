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
seq_n <- 300#c(25,50,100,200,300,500)
seq_p <- c(5,10,15,20)
seq_lambda <- c(1,5,10,15)#seq(0,50,5)#
seq_rep <- 1#1:5
Allbeta <- beta <- c(1:5,rep(0, max(seq_p) - 5))# rbinom(p, size = 1, prob = 1/seq(1,p/2, by = 0.5))
iter.max <- 5000

grid <- expand.grid(n =  seq_n,
                    p =  seq_p,
                    lambda =  seq_lambda,
                    rep =  seq_rep)
n.grid <- nrow(grid)

export.fit <- c("cv","iteration","maxDiff","minDiff","time")
names.results <- c("n.obs","n.param","lambda","rep",
                   "lvm.time","L1.zeroBeta",
                   #paste("L1",export.fit, sep = "."),
                   paste("L1ISTA",export.fit, sep = "."),
                   paste("L1FISTA",export.fit, sep = "."),
                   paste("L2",export.fit, sep = "."),
                   paste("L2ISTA",export.fit, sep = "."),
                   paste("L2FISTA",export.fit, sep = ".")
)
dt.results <- data.table(matrix(as.numeric(NA), nrow = n.grid, ncol = length(names.results)))
names(dt.results) <- names.results

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
  plvm.model <- penalize(lvm.model)
  
  ## lava
  tps1 <- system.time(
  elvm.fit <- estimate(lvm.model, df.data, 
                       control = list(constrain = TRUE))
  )

#   tps2 <- system.time(
#   elvm1.fit <- estimate(plvm.model, data = df.data, lambda1 = lambda1, 
#                         control = list(constrain = TRUE, iter.max = iter.max))
#   )
  tps2bis <- system.time(
    elvm1ISTA.fit <- estimate(plvm.model,  data = df.data, lambda1 = lambda1, method = "ISTA",
                              control = list(constrain = TRUE, iter.max = iter.max))
  )
  tps2ter <- system.time(
    elvm1FISTA.fit <- estimate(plvm.model,  data = df.data, lambda1 = lambda1,  method = "FISTA",
                               control = list(constrain = TRUE, iter.max = iter.max))
  )
  tps3 <- system.time(
    elvm2.fit <- estimate(plvm.model,  data = df.data, lambda2 = lambda2, method = "NR",
                          control = list(constrain = TRUE, iter.max = iter.max))
  )
#   tps3bis <- system.time(
#     elvm2ISTA.fit <- estimate(plvm.model,  data = df.data, lambda2 = lambda2, method = "ISTA",
#                               control = list(constrain = TRUE, iter.max = iter.max))
#   )
#   tps3ter <- system.time(
#     elvm2FISTA.fit <- estimate(plvm.model,  data = df.data, lambda2 = lambda2,  method = "FISTA",
#                                control = list(constrain = TRUE, iter.max = iter.max))
#   )
  
    
  ## penalized
  penalized1.fit <- penalized:::penalized(response = df.data$Y, penalized = df.data[,all.vars(formula_tempo)[-1]], 
                                          lambda1 = lambda1 * coef(elvm1FISTA.fit)["Y,Y"], lambda2 = 0, trace = FALSE )
  
  penalized2.fit <- penalized:::penalized(response = df.data$Y, penalized = df.data[,all.vars(formula_tempo)[-1]], 
                                          lambda1 = 0, lambda2 = lambda2 * coef(elvm2.fit)["Y,Y"], trace = FALSE )
  
  
  dt.results[iter_grid, c("n.obs","n.param","lambda","rep") := .(n,p,grid[iter_grid,"lambda"],grid[iter_grid,"rep"])]
  dt.results[iter_grid, c("lvm.time","L1ISTA.time","L1FISTA.time","L2.time") := .(tps1[3],tps2bis[3],tps2ter[3],tps3[3])]
  
  param_penalized1 <- c(penalized1.fit@unpenalized, penalized1.fit@penalized, NA) #penalized1.fit@nuisance$sigma2)
  dt.results[iter_grid, L1.zeroBeta := sum( abs(penalized1.fit@penalized) < 1e-10)]
#   dt.results[iter_grid, L1.cv := if("convergence" %in% names(elvm1.fit$opt)){elvm1.fit$opt$convergence}else{NA} ]
#   dt.results[iter_grid, L1.iteration := if("iterations" %in% names(elvm1.fit$opt)){elvm1.fit$opt$iterations}else{NA} ]
#   dt.results[iter_grid, L1.maxDiff := max(abs(param_penalized1-coef(elvm1.fit)[1:(p+1)]))]
#   dt.results[iter_grid, L1.minDiff := min(abs(param_penalized1-coef(elvm1.fit)[1:(p+1)]))]
  
  dt.results[iter_grid, L1ISTA.cv := if("convergence" %in% names(elvm1ISTA.fit$opt)){elvm1ISTA.fit$opt$convergence}else{NA} ]
  dt.results[iter_grid, L1ISTA.iteration := if("iterations" %in% names(elvm1ISTA.fit$opt)){elvm1ISTA.fit$opt$iterations}else{NA} ]
  dt.results[iter_grid, L1ISTA.maxDiff := max(abs(param_penalized1-coef(elvm1ISTA.fit)), na.rm = TRUE)]
  dt.results[iter_grid, L1ISTA.minDiff := min(abs(param_penalized1-coef(elvm1ISTA.fit)), na.rm = TRUE)]
  
  dt.results[iter_grid, L1FISTA.cv := if("convergence" %in% names(elvm1FISTA.fit$opt)){elvm1FISTA.fit$opt$convergence}else{NA} ]
  dt.results[iter_grid, L1FISTA.iteration := if("iterations" %in% names(elvm1FISTA.fit$opt)){elvm1FISTA.fit$opt$iterations}else{NA} ]
  dt.results[iter_grid, L1FISTA.maxDiff := max(abs(param_penalized1-coef(elvm1FISTA.fit)), na.rm = TRUE)]
  dt.results[iter_grid, L1FISTA.minDiff := min(abs(param_penalized1-coef(elvm1FISTA.fit)), na.rm = TRUE)]
  
  param_penalized2 <- c(penalized2.fit@unpenalized, penalized2.fit@penalized, NA) #penalized2.fit@nuisance$sigma2)
  dt.results[iter_grid, L2.cv := if("convergence" %in% names(elvm2.fit$opt)){elvm2.fit$opt$convergence}else{NA} ]
  dt.results[iter_grid, L2.iteration := if("iterations" %in% names(elvm2.fit$opt)){elvm2.fit$opt$iterations}else{NA} ]
  dt.results[iter_grid, L2.maxDiff := max(abs(param_penalized2-coef(elvm2.fit)), na.rm = TRUE)]
  dt.results[iter_grid, L2.minDiff := min(abs(param_penalized2-coef(elvm2.fit)), na.rm = TRUE)]
  
#   dt.results[iter_grid, L2ISTA.cv := if("convergence" %in% names(elvm2ISTA.fit$opt)){elvm2ISTA.fit$opt$convergence}else{NA} ]
#   dt.results[iter_grid, L2ISTA.iteration := if("iterations" %in% names(elvm2ISTA.fit$opt)){elvm2ISTA.fit$opt$iterations}else{NA} ]
#   dt.results[iter_grid, L2ISTA.maxDiff := max(abs(param_penalized2-coef(elvm2ISTA.fit)), na.rm = TRUE)]
#   dt.results[iter_grid, L2ISTA.minDiff := min(abs(param_penalized2-coef(elvm2ISTA.fit)), na.rm = TRUE)]
#   
#   dt.results[iter_grid, L2FISTA.cv := if("convergence" %in% names(elvm2FISTA.fit$opt)){elvm2FISTA.fit$opt$convergence}else{NA} ]
#   dt.results[iter_grid, L2FISTA.iteration := if("iterations" %in% names(elvm2FISTA.fit$opt)){elvm2FISTA.fit$opt$iterations}else{NA} ]
#   dt.results[iter_grid, L2FISTA.maxDiff := max(abs(param_penalized2-coef(elvm2FISTA.fit)), na.rm = TRUE)]
#   dt.results[iter_grid, L2FISTA.minDiff := min(abs(param_penalized2-coef(elvm2FISTA.fit)), na.rm = TRUE)]
  
}
print(dt.results)

# save(dt.results, file = paste0("C:/Users/hpl802/Documents/Projects/LVM/CV-test/dt_results.RData"))


# dt.results[, lapply(.SD,mean), .SDcols = c("L1ISTA.cv","L1FISTA.cv","L2.cv","L2ISTA.cv","L2FISTA.cv")]
# dt.results[, lapply(.SD,quantile, na.rm = TRUE), .SDcols = c("L1ISTA.maxDiff","L1FISTA.maxDiff","L2.maxDiff","L2ISTA.maxDiff","L2FISTA.maxDiff")]
# dt.results[, lapply(.SD,quantile), .SDcols = c("L1ISTA.iteration","L1FISTA.iteration","L2.iteration","L2ISTA.iteration","L2FISTA.iteration")]
# dt.results[, lapply(.SD,quantile, na.rm = TRUE), .SDcols = c("L1ISTA.time","L1FISTA.time","L2.time","L2ISTA.time","L2FISTA.time")]

# dt.results[ lambda == 1, lapply(.SD,quantile), .SDcols = c("L1ISTA.maxDiff","L1FISTA.maxDiff","L2.maxDiff","L2ISTA.maxDiff","L2FISTA.maxDiff")]

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

