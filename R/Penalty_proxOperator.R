proxL1 <- function(x, step, lambda, test.penalty){
 if(test.penalty){
    max(0, (abs(x) - lambda * step)) * sign(x)
  }else{
    x
  }
}

proxL2 <- function(x, step, lambda, test.penalty){
  if(test.penalty){
    1 / (1 + lambda * step) * x
  }else{
    x
  }
}

proxE2 <- function(x, step, lambda){
  max(0, 1- sqrt(length(x)) * lambda * step/norm(x, type = "2"))*x
}


#### need to adapt test.penalty
# x <- 1:10
# step <- 0.05
# lambda <- 10
# test.penalty <- abs(round(rnorm(10, sd = 2)))
# proxE2(x, step, lambda, test.penalty)

