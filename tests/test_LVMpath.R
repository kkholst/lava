set.seed(10)
n <- 500
lvm.model <- lvm(list(y~X1+X2+X3+X4))
df.data <- sim(lvm.model,n)
plvm.model <- penalize(lvm.model)

lvm.fit <- estimate(lvm.model,  data = df.data,
                      control = list(constrain = TRUE))

#### look at sigma
lambda.lvm <- 0.5
p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = lambda.lvm,
                      control = list(constrain = TRUE, iter.max = 1000))
p1lvm.fit
coef(p1lvm.fit)
lambda.penalized <- lambda.lvm*n*coef(p1lvm.fit)["y,y"]
penalized.fit <- penalized(df.data$y,df.data[,c("X1", "X2", "X3", "X4")], 
                 lambda1 = lambda.penalized)
coef(p1lvm.fit) - c(penalized.fit@unpenalized, penalized.fit@penalized,penalized.fit@nuisance$sigma2)
lambda.lvm-lambda.penalized/(n*penalized.fit@nuisance$sigma2)


#### check the whole path
path.fit <- penalized(df.data$y, df.data[,c("X1", "X2", "X3", "X4")], 
                      steps = "Park" )
seq_lambda.penalized <- unlist(lapply(path.fit, function(x){x@lambda1}))
seq_sigma.penalized <- unlist(lapply(path.fit, function(x){x@nuisance$sigma2}))
seq_lambda.lvm <- seq_lambda.penalized / (n * seq_sigma.penalized)

iter_path <- 5
p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = seq_lambda.lvm[iter_path], regularizationPath = TRUE,
                      control = list(constrain = TRUE, iter.max = 1000))
coef(p1lvm.fit)
c(path.fit[[iter_path]]@unpenalized, path.fit[[iter_path]]@penalized, path.fit[[iter_path]]@nuisance$sigma2)

penalized.fit <- penalized(df.data$y,df.data[,c("X1", "X2", "X3", "X4")], 
                           lambda1 = 356.24646129)

coef(p1lvm.fit <- estimate(plvm.model,  data = df.data, lambda1 = 117 / n, regularizationPath = FALSE,
                      control = list(constrain = TRUE, iter.max = 1000)))
99.1434369149134 119.082673005792 500.032469917492  

coef(p1lvm.fit)

coef(res[[iter_path]])


seq_lambda



lambda1 <- rep(0, length(start))
lambda1[penalty$index.coef] <- 2000
do.call(control$proxGrad.method,
               list(start = start, step = step, proxOperator = proxOperator, gradient = gradient, 
                    lambda1 = lambda1, lambda2 = lambda2,
                    iter.max = control$iter.max, abs.tol = control$abs.tol, rel.tol = control$rel.tol))$par

res$objective <- objective(res$par, penalty = penalty)
