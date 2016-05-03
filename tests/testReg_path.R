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

validLVM <- function(x){
  library(penalized)
  
  name.outcome <- strsplit(x$penalty$names.varCoef, split = ",", fixed = TRUE)[[1]][1]
  resPenalized <- penalized(as.formula(paste0(name.outcome,"~.")), lambda1 = x$penalty$lambda1*coef(x)[x$penalty$names.varCoef], data = x$data$model.frame)
  diffCoef <- coef(x) - coef2.penalized(resPenalized)
  
  return(diffCoef)
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
  
#   elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
#                              fix.sigma = TRUE, lambda1 = penalized.PathL1[[iter_l]]@lambda1,
#                              control = list(constrain = TRUE, iter.max = 1000, fast = 0, trace = TRUE))
#   print( range(coef(elvm.fit_tempo) - coef2.penalized(penalized.PathL1, iter_lambda = iter_l)) )
#   print( elvm.fit_tempo$opt$iterations )
#   elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
#                              lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
#                              control = list(constrain = TRUE, iter.max = 1000, fast = 0, trace = FALSE))
#   print( range(coef(elvm.fit_tempo) - coef2.penalized(penalized.PathL1, iter_lambda = iter_l)) )
#   print( elvm.fit_tempo$opt$iterations )
  
  elvm.fit_tempo <- estimate(plvm.model,  data = df.data, 
                             lambda1 = penalized.PathL1[[iter_l]]@lambda1/penalized.PathL1[[iter_l]]@nuisance$sigma2,
                             control = list(constrain = TRUE, iter.max = 1000, fast = 3, trace = FALSE))
  print( range(coef(elvm.fit_tempo) - coef2.penalized(penalized.PathL1, iter_lambda = iter_l)) )
  print( elvm.fit_tempo$opt$iterations )
}

#### check path ####
elvm.PathL1_fixed <- estimate(plvm.model,  data = df.data, 
                              fix.sigma = TRUE, regularizationPath = TRUE,
                              control = list(constrain = TRUE, iter.max = 5000))

elvm.PathL1_fixed$opt$message
coef2.penalized(penalized.PathL1)

beta <- coef(estimate(lvm.model, df.data))
elvm.PathL1_free <- estimate(plvm.model,  data =  df.data, 
                              regularizationPath = TRUE, lambda2 = 0,
                              control = list(constrain = FALSE, iter.max = 5000, start = coef(estimate(lvm.model, df.data))))
elvm.PathL1_free$opt$message

# list(elvm.PathL1_fixed$opt$message$lambda1,
#      elvm.PathL1_free$opt$message$lambda1 * elvm.PathL1_free$opt$message[,"Y,Y"])
# [1] 1617.81111 1617.64933 1002.92748  606.39734  606.33670   29.44781   29.44487   24.99796    0.00000

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

#### 3- Group Lasso ####
# gglasso
# grplasso
# grpreg

####
library(gglasso)
data(bardet)
group1 <- rep(1:4,each=1)

bardet$x <- bardet$x[,1:5]

# fit group lasso penalized least squares
m1 <- gglasso(x=bardet$x,y=bardet$y,group=group1,loss="ls")
m1$b0
m1$beta
m1$npasses

df.bardet <- data.frame(bardet)
names(df.bardet) <- gsub(".","",names(df.bardet), fixed = TRUE)
lvm.model2 <- lvm(as.formula(paste0("y ~ ", paste(names(df.bardet)[1:ncol(df.bardet)], collapse = "+"))))
          
plvm.model2 <- penalize(lvm.model2) 
plvm.model2$penalty$group.penaltyCoef <- group1

eplvm.model2 <- estimate(plvm.model2,  data = df.bardet,
                         lambda1 = 75,
                         control = list(constrain = FALSE, iter.max = 1000, step = NULL, trace = TRUE))
coef(eplvm.model2)
####
library(grplasso)
data(splice)

## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

## Fit a logistic model 
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 20,
                       contrasts = contr, center = TRUE, standardize = TRUE)

####
library(grpreg)
data(birthwt.grpreg)
X <- as.matrix(birthwt.grpreg[,-1:-2])
y <- birthwt.grpreg$bwt
group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
fit <- grpreg(X,y,group,penalty="grLasso")
plot(fit)

#### 4- Several regressions #####
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
