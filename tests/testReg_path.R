rm(list = ls())

library(penalized)
library(optimx)
library(numDeriv)
library(data.table)

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

#### 1- only L1 ####
set.seed(10)
n <- 500
formula.lvm <- as.formula(paste0("Y~",paste(paste0("X",1:5), collapse = "+")))
lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm) <- as.list( c(rep(0,2),1:3) )
distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(formula.lvm)
plvm.model <- penalize(lvm.model)


####
penalized.PathL1 <- penalized(Y ~  ., data = df.data, steps = "Park", trace = TRUE)

seq_lambda <- unlist(lapply(penalized.PathL1, function(x){x@lambda1}))

#### check fix lambda ####
for(iter_l in 1:length(seq_lambda)){
  cat(iter_l, " : ",penalized.PathL1[[iter_l]]@lambda1,"\n")
  
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
                             fix.sigma = TRUE, lambda1 = penalized.PathL1[[iter_l]]@lambda1,
                             control = list(constrain = TRUE, iter.max = 5000))
  print( range(coef(elvm.fit_tempo) - coef2.penalized(penalized.PathL1, iter_lambda = iter_l)) )
  
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                             control = list(constrain = TRUE, iter.max = 5000))
  print( range(coef(elvm.fit_tempo) - coef2.penalized(penalized.PathL1, iter_lambda = iter_l)) )
}

#### check path ####
elvm.PathL1_fixed <- estimate(plvm.model,  data = df.data, 
                              fix.sigma = TRUE, regularizationPath = TRUE,
                              control = list(constrain = TRUE, iter.max = 5000))

elvm.PathL1_fixed$opt$message
coef2.penalized(penalized.PathL1)


elvm.PathL1_free <- estimate(plvm.model,  data =  df.data, 
                              fix.sigma = FALSE, regularizationPath = TRUE, lambda2 = 0,
                              control = list(constrain = FALSE, iter.max = 5000))
elvm.PathL1_free$opt$message

# list(elvm.PathL1_fixed$opt$message$lambda1,
#      elvm.PathL1_free$opt$message$lambda1 * elvm.PathL1_free$opt$message[,"Y,Y"])

iter_path <- 5
test <- estimate(plvm.model,  data = df.data, 
                 fix.sigma = FALSE, lambda1 = elvm.PathL1_free$opt$message[iter_path,"lambda1"], lambda2 = 0,#,
                 control = list(constrain = TRUE, iter.max = 5000))
rbind(coef(test),
      elvm.PathL1_free$opt$message[iter_path,-(1:2)])

#### 2- L1 and L2 ####
penalized.PathL12 <- penalized(Y ~  ., data = df.data, steps = "Park", lambda2 = 100, trace = FALSE)

elvm.PathL12_fixed <- estimate(plvm.model,  data = df.data, 
                              fix.sigma = TRUE, regularizationPath = TRUE, lambda2 = 100,
                              control = list(constrain = FALSE, iter.max = 5000))
elvm.PathL12_fixed$opt$message
coef2.penalized(penalized.PathL12)


elvm.PathL12_free <- estimate(plvm.model,  data = df.data, 
                             fix.sigma = FALSE, regularizationPath = TRUE, lambda2 = 100,
                             control = list(constrain = FALSE, iter.max = 5000))
elvm.PathL12_free$opt$message

iter_path <- 5
test <- estimate(plvm.model,  data = df.data, 
                 fix.sigma = FALSE, lambda1 = elvm.PathL12_free$opt$message[iter_path,"lambda1"], lambda2 = 100,#,
                 control = list(constrain = TRUE, iter.max = 5000))
rbind(coef(test),
      elvm.PathL12_free$opt$message[iter_path,-(1:2)])

#### 3- Several regressions #####
set.seed(10)
n <- 500
formula.lvm1 <- as.formula(paste0("Y1~",paste(paste0("X",1:5), collapse = "+")))
formula.lvm2 <- as.formula(paste0("Y2~",paste(paste0("X",4:6), collapse = "+")))
#formula.lvm3 <- as.formula(paste0("Y3~",paste(paste0("X",5:10), collapse = "+")))

lvm.modelSim <- lvm()
regression(lvm.modelSim, formula.lvm1) <- as.list( c(rep(0,2),1:3) )
regression(lvm.modelSim, formula.lvm2) <- as.list( c(rep(0,1),1:2) )
#regression(lvm.modelSim, formula.lvm3) <- as.list( c(rep(0,4),2:3) )

# distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
# distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
# distribution(lvm.modelSim, ~Y) <- normal.lvm(sd = 2)
df.data <- sim(lvm.modelSim,n)

lvm.model <- lvm(list(formula.lvm1, formula.lvm2))
plvm.model <- penalize(lvm.model)


penalized.PathL1 <- penalized(formula.lvm1, data = df.data, steps = "Park", trace = TRUE)
coef2.penalized(penalized.PathL1)

elvm.PathL1_fixed <- estimate(plvm.model,  data = df.data, 
                              fix.sigma = FALSE, regularizationPath = TRUE,
                              control = list(constrain = FALSE, iter.max = 5000))
Mres <- elvm.PathL1_fixed$opt$message
Mres$lambda1
Mres[,c(1:2,grep("Y1",names(Mres),fixed = TRUE))]


lvm.test <- estimate(plvm.model,  data = df.data, 
                     fix.sigma = FALSE, lambda1 = Mres$lambda1[10],
                     control = list(constrain = FALSE, iter.max = 5000))
coef(lvm.test)
