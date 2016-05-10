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

proxE2 <- function(x, step, lambda){ # adapted from Simon 2013
  max(0, 1- sqrt(length(x)) * lambda * step/norm(x, type = "2"))*x
}

