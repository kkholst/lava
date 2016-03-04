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


