# Calculate simultaneous confidence limits by Scheffe's method

Function to compute the Scheffe corrected confidence interval for the
regression line

## Usage

``` r
scheffe(model, newdata = model.frame(model), level = 0.95)
```

## Arguments

- model:

  Linear model

- newdata:

  new data frame

- level:

  confidence level (0.95)

## Examples

``` r
x <- rnorm(100)
d <- data.frame(y=rnorm(length(x),x),x=x)
l <- lm(y~x,d)
plot(y~x,d)
abline(l)
d0 <- data.frame(x=seq(-5,5,length.out=100))
d1 <- cbind(d0,predict(l,newdata=d0,interval="confidence"))
d2 <- cbind(d0,scheffe(l,d0))
lines(lwr~x,d1,lty=2,col="red")
lines(upr~x,d1,lty=2,col="red")
lines(lwr~x,d2,lty=2,col="blue")
lines(upr~x,d2,lty=2,col="blue")
```
