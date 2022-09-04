
# Table of Contents

1.  [Estimation](#org2fd4bce)
2.  [Cross-validation](#orgac174ab)
3.  [Specification of general functional forms](#org8a61fb4)
4.  [Mixture models](#orgedc0730)

<!-- nonlinear.Rmd is generated from nonlinear.org. Please edit that file -->

    knitr::opts_chunk$set(
      collapse = TRUE,
      comment = "#>"
    )
    mets <- lava:::versioncheck('mets', 1)
    fullVignette <- Sys.getenv("_R_FULL_VIGNETTE_") %in% c("1","TRUE") ||
      length(list.files("data"))==0

    library('lava')

We consider the measurement models given by

$$X_{j} = \eta_{1} + \epsilon_{j}^{x}, \quad j=1,2,3$$
$$Y_{j} = \eta_{2} + \epsilon_{j}^{y}, \quad j=1,2,3$$
and with a structural model given by
$$\eta_{2} = f(\eta_{1}) + Z + \zeta_{2}\label{ex:eta2}$$
$$\eta_{1} = Z + \zeta_{1}\label{ex:eta1}$$
with iid measurement errors
$\epsilon_{j}^{x},\epsilon_{j}^{y},\zeta_{1},\zeta_{2}\sim\mathcal{N}(0,1),
j=1,2,3.$ and standard normal distributed covariate $Z$.  To
simulate from this model we use the following syntax:

    f <- function(x) cos(1.25*x) + x - 0.25*x^2
    m <- lvm(x1+x2+x3 ~ eta1, y1+y2+y3 ~ eta2, latent=~eta1+eta2)
    regression(m) <- eta1+eta2 ~ z
    functional(m, eta2~eta1) <- f
    
    d <- sim(m, n=200, seed=42) # Default is all parameters are 1

    plot(m, plot.engine="visNetwork")

We refer to <holst_budtzjorgensen_2013> for details on the syntax for
model specification.


<a id="org2fd4bce"></a>

# Estimation

To estimate the parameters using the two-stage estimator described in <holst_budtzjorgensen_2020>,  the first step is now to specify the measurement models

    m1 <- lvm(x1+x2+x3 ~ eta1, eta1 ~ z, latent=~eta1)
    m2 <- lvm(y1+y2+y3 ~ eta2, eta2 ~ z, latent=~eta2)

Next, we specify a quadratic relationship between the two latent variables

    nonlinear(m2, type="quadratic") <- eta2 ~ eta1

and the model can then be estimated using the two-stage estimator

    e1 <- twostage(m1, m2, data=d)
    e1

                        Estimate Std. Error  Z-value  P-value
    Measurements:                                            
       y2~eta2           0.97686    0.08406 11.62146   <1e-12
       y3~eta2           1.04485    0.04531 23.06093   <1e-12
    Regressions:                                             
       eta2~z            0.88513   16.11882  0.05491   0.9562
       eta2~eta1_1       1.14072   18.98292  0.06009   0.9521
       eta2~eta1_2      -0.45055    9.16379 -0.04917   0.9608
    Intercepts:                                              
       y2               -0.12198    0.10935 -1.11549   0.2646
       y3               -0.09879    0.10569 -0.93471   0.3499
       eta2              0.67814   15.00825  0.04518    0.964
    Residual Variances:                                      
       y1                1.30730    0.48572  2.69149         
       y2                1.11056    1.07092  1.03702         
       y3                0.80961    0.76919  1.05255         
       eta2              2.08483    7.57087  0.27537

We see a clear statistically significant effect of the second order
term (`eta2~eta1_2`). For comparison we can also estimate the full MLE
of the linear model:

    e0 <- estimate(regression(m1%++%m2, eta2~eta1), d)
    estimate(e0,keep="^eta2~[a-z]",regex=TRUE) ## Extract coef. matching reg.ex.

              Estimate Std.Err    2.5% 97.5%   P-value
    eta2~eta1   1.4140  0.2261 0.97083 1.857 4.014e-10
    eta2~z      0.6374  0.2778 0.09291 1.182 2.177e-02

Next, we calculate predictions from the quadratic model using the estimated parameter coefficients
$$
\mathbb{E}_{\widehat{\theta}_{2}}(\eta_{2} \mid \eta_{1}, Z=0),
$$

    newd <- expand.grid(eta1=seq(-4, 4, by=0.1), z=0)
    pred1 <- predict(e1, newdata=newd, x=TRUE)
    head(pred1)

                 y1         y2         y3       eta2
    [1,] -11.093569 -10.958869 -11.689950 -11.093569
    [2,] -10.623561 -10.499736 -11.198861 -10.623561
    [3,] -10.162565 -10.049406 -10.717187 -10.162565
    [4,]  -9.710579  -9.607878 -10.244928  -9.710579
    [5,]  -9.267605  -9.175153  -9.782084  -9.267605
    [6,]  -8.833641  -8.751230  -9.328656  -8.833641

To obtain a potential better fit we next proceed with a natural cubic spline

    kn <- seq(-3,3,length.out=5)
    nonlinear(m2, type="spline", knots=kn) <- eta2 ~ eta1
    e2 <- twostage(m1, m2, data=d)
    e2

                        Estimate Std. Error  Z-value  P-value
    Measurements:                                            
       y2~eta2           0.97752    0.09601 10.18187   <1e-12
       y3~eta2           1.04508    0.08418 12.41419   <1e-12
    Regressions:                                             
       eta2~z            0.86729   15.75435  0.05505   0.9561
       eta2~eta1_1       2.86231   21.41242  0.13368   0.8937
       eta2~eta1_2       0.00344    3.46229  0.00099   0.9992
       eta2~eta1_3      -0.26270   11.47890 -0.02289   0.9817
       eta2~eta1_4       0.50778   14.09762  0.03602   0.9713
    Intercepts:                                              
       y2               -0.12185    0.10980 -1.10974   0.2671
       y3               -0.09874    0.10627 -0.92918   0.3528
       eta2              1.83814   69.73223  0.02636    0.979
    Residual Variances:                                      
       y1                1.31286    0.49428  2.65613         
       y2                1.10412    1.03272  1.06914         
       y3                0.81124    0.72916  1.11257         
       eta2              1.99404    6.88395  0.28966

Confidence limits can be obtained via the Delta method using the `estimate` method:

    p <- cbind(eta1=newd$eta1,
    	  estimate(e2,f=function(p) predict(e2,p=p,newdata=newd))$coefmat)
    head(p)

       eta1  Estimate   Std.Err      2.5%    97.5%   P-value
    p1 -4.0 -9.611119 106.92612 -219.1825 199.9602 0.9283781
    p2 -3.9 -9.324887 105.30963 -215.7280 197.0782 0.9294417
    p3 -3.8 -9.038656 103.71215 -212.3107 194.2334 0.9305512
    p4 -3.7 -8.752425 102.13458 -208.9325 191.4277 0.9317089
    p5 -3.6 -8.466193 100.57786 -205.5952 188.6628 0.9329169
    p6 -3.5 -8.179962  99.04296 -202.3006 185.9407 0.9341775

The fitted function can be obtained with the following code:

    plot(I(eta2-z) ~ eta1, data=d, col=Col("black",0.5), pch=16,
         xlab=expression(eta[1]), ylab=expression(eta[2]), xlim=c(-4,4))
    lines(Estimate ~ eta1, data=as.data.frame(p), col="darkblue", lwd=5)
    confband(p[,1], lower=p[,4], upper=p[,5], polygon=TRUE,
    	 border=NA, col=Col("darkblue",0.2))


<a id="orgac174ab"></a>

# Cross-validation

A more formal comparison of the different models can be obtained by
cross-validation. Here we specify linear, quadratic and cubic spline
models with 4 and 9 degrees of freedom.

    m2a <- nonlinear(m2, type="linear", eta2~eta1)
    m2b <- nonlinear(m2, type="quadratic", eta2~eta1)
    kn1 <- seq(-3,3,length.out=5)
    kn2 <- seq(-3,3,length.out=8)
    m2c <- nonlinear(m2, type="spline", knots=kn1, eta2~eta1)
    m2d <- nonlinear(m2, type="spline", knots=kn2, eta2~eta1)

To assess the model fit average RMSE is estimated with 5-fold
cross-validation repeated two times

    ## Scale models in stage 2 to allow for a fair RMSE comparison
    d0 <- d
    for (i in endogenous(m2))
        d0[,i] <- scale(d0[,i],center=TRUE,scale=TRUE)
    ## Repeated 5-fold cross-validation:
    ff <- lapply(list(linear=m2a,quadratic=m2b,spline4=m2c,spline6=m2d),
    	    function(m) function(data,...) twostage(m1,m,data=data,stderr=FALSE,control=list(start=coef(e0),contrain=TRUE)))
    fit.cv <- targeted::cv(ff,data=d,K=5,rep=2,mc.cores=parallel::detectCores(),seed=1)

    Error in `[.data.frame`(response, fold[[k]]) : undefined columns selected

    ## To save time building the vignettes on CRAN, we cache time consuming computations
    if (fullVignette) {
      fit.cv$fit <- NULL
      saveRDS(fit.cv, "data/nonlinear_fitcv.rds", version=2)
    } else {
      fit.cv <- readRDS("data/nonlinear_fitcv.rds")
    }

    Error in fit.cv$fit <- NULL : object 'fit.cv' not found

    summary(fit.cv)

    Error in summary(fit.cv) : object 'fit.cv' not found

Here the RMSE is in favour of the splines model with 4 degrees of freedom:

    fit <- lapply(list(m2a,m2b,m2c,m2d),
    	     function(x) {
    		 e <- twostage(m1,x,data=d)
    		 pr <- cbind(eta1=newd$eta1,predict(e,newdata=newd$eta1,x=TRUE))
    		 return(list(estimate=e,predict=as.data.frame(pr)))
    	     })
    
    plot(I(eta2-z) ~ eta1, data=d, col=Col("black",0.5), pch=16,
         xlab=expression(eta[1]), ylab=expression(eta[2]), xlim=c(-4,4))
    col <- c("orange","darkred","darkgreen","darkblue")
    lty <- c(3,4,1,5)
    for (i in seq_along(fit)) {
        with(fit[[i]]$pr, lines(eta2 ~ eta1, col=col[i], lwd=4, lty=lty[i]))
    }
    legend("bottomright",
      c("linear","quadratic","spline(df=4)","spline(df=6)"),
      col=col, lty=lty, lwd=3)

For convenience, the function `twostageCV` can be used to do the
cross-validation (also for choosing the mixture distribution via the \`\`nmix\`\` argument, see the section
below). For example,

    selmod <- twostageCV(m1, m2, data=d, df=2:4, nmix=1:2,
    	    nfolds=2, rep=1, mc.cores=parallel::detectCores())

    ## To save time building the vignettes on CRAN, we cache time consuming computations
    if (fullVignette) {
      saveRDS(summary(selmod), "data/nonlinear_selmod.rds", version=2)
    } else {
      selmod <- readRDS("data/nonlinear_selmod.rds")
    }

    Error in summary(selmod) : object 'selmod' not found

applies cross-validation (here just 2 folds for simplicity) to select the best splines with
degrees of freedom varying from from 1-3 (the linear model is
automatically included)

    selmod

    Error: object 'selmod' not found


<a id="org8a61fb4"></a>

# Specification of general functional forms

Next, we show how to specify a general functional relation of
multiple different latent or exogenous variables. This is achieved via
the `predict.fun` argument. To illustrate this we include interactions
between the latent variable $\eta_{1}$ and a dichotomized version of
the covariate $z$

    d$g <- (d$z<0)*1 ## Group variable
    mm1 <- regression(m1, ~g)  # Add grouping variable as exogenous variable (effect specified via 'predict.fun')
    mm2 <- regression(m2, eta2~ u1+u2+u1:g+u2:g+z)
    pred <- function(mu,var,data,...) {
        cbind("u1"=mu[,1],"u2"=mu[,1]^2+var[1],
    	  "u1:g"=mu[,1]*data[,"g"],"u2:g"=(mu[,1]^2+var[1])*data[,"g"])
    }
    ee1 <- twostage(mm1, model2=mm2, data=d, predict.fun=pred)
    estimate(ee1,keep="eta2~u",regex=TRUE)

              Estimate Std.Err   2.5% 97.5% P-value
    eta2~u1     0.9891  17.805 -33.91 35.89  0.9557
    eta2~u2    -0.3962   9.971 -19.94 19.15  0.9683
    eta2~u1:g   0.4487  13.679 -26.36 27.26  0.9738
    eta2~u2:g   0.0441   7.178 -14.02 14.11  0.9951

A formal test show no statistically significant effect of this interaction

    summary(estimate(ee1,keep="(:g)", regex=TRUE))

    Call: estimate.default(x = ee1, keep = "(:g)", regex = TRUE)
    __________________________________________________
              Estimate Std.Err   2.5% 97.5% P-value
    eta2~u1:g   0.4487  13.679 -26.36 27.26  0.9738
    eta2~u2:g   0.0441   7.178 -14.02 14.11  0.9951
    
     Null Hypothesis: 
      [eta2~u1:g] = 0
      [eta2~u2:g] = 0 
     
    chisq = 0.0016, df = 2, p-value = 0.9992


<a id="orgedc0730"></a>

# Mixture models

Lastly, we demonstrate how the distributional assumptions of stage 1
model can be relaxed by letting the conditional distribution of the
latent variable given covariates follow a Gaussian mixture
distribution. The following code explictly defines the parameter
constraints of the model by setting the intercept of the first
indicator variable, $x_{1}$, to zero and the factor loading
parameter of the same variable to one.

    m1 <- baptize(m1)  ## Label all parameters
    intercept(m1, ~x1+eta1) <- list(0,NA) ## Set intercept of x1 to zero. Remove the label of eta1
    regression(m1,x1~eta1) <- 1 ## Factor loading fixed to 1

The mixture model may then be estimated using the `mixture` method
(note, this requires the `mets` package to be installed), where the
Parameter names shared across the different mixture components given
in the `list` will be constrained to be identical in the mixture
model. Thus, only the intercept of $\eta_{1}$ is allowed to vary
between the mixtures.

    set.seed(1)
    em0 <- mixture(m1, k=2, data=d)

To decrease the risk of using a local maximizer of the likelihood we
can rerun the estimation with different random starting values

    em0 <- NULL
    ll <- c()
    for (i in 1:5) {
        set.seed(i)
        em <- mixture(m1, k=2, data=d, control=list(trace=0))
        ll <- c(ll,logLik(em))
        if (is.null(em0) || logLik(em0)<tail(ll,1))
    	em0 <- em
    }

    ## To save time building the vignettes on CRAN, we cache time consuming computations
    if (fullVignette) {
      saveRDS(em0, "data/nonlinear_em0.rds", version=2)
    } else {
      em0 <- readRDS("data/nonlinear_em0.rds")
    }

    Error in gzfile(file, mode) : cannot open the connection

    summary(em0)

    Cluster 1 (n=162, Prior=0.776):
    --------------------------------------------------
                        Estimate Std. Error Z value  Pr(>|z|)
    Measurements:                                            
       x1~eta1           1.00000                             
       x2~eta1           0.99581  0.07940   12.54098   <1e-12
       x3~eta1           1.06344  0.08436   12.60538   <1e-12
    Regressions:                                             
       eta1~z            1.06675  0.08527   12.50987   <1e-12
    Intercepts:                                              
       x1                0.00000                             
       x2                0.03845  0.09890    0.38883 0.6974  
       x3               -0.02549  0.10333   -0.24667 0.8052  
       eta1              0.20923  0.13162    1.58961 0.1119  
    Residual Variances:                                      
       x1                0.98540  0.13316    7.40021         
       x2                0.97181  0.13156    7.38694         
       x3                1.01316  0.14294    7.08810         
       eta1              0.29048  0.11129    2.61001         
    
    Cluster 2 (n=38, Prior=0.224):
    --------------------------------------------------
                        Estimate Std. Error Z value  Pr(>|z|) 
    Measurements:                                             
       x1~eta1           1.00000                              
       x2~eta1           0.99581  0.07940   12.54098   <1e-12 
       x3~eta1           1.06344  0.08436   12.60538   <1e-12 
    Regressions:                                              
       eta1~z            1.06675  0.08527   12.50987   <1e-12 
    Intercepts:                                               
       x1                0.00000                              
       x2                0.03845  0.09890    0.38883 0.6974   
       x3               -0.02549  0.10333   -0.24667 0.8052   
       eta1             -1.44295  0.25868   -5.57806 2.432e-08
    Residual Variances:                                       
       x1                0.98540  0.13316    7.40021          
       x2                0.97181  0.13156    7.38694          
       x3                1.01316  0.14294    7.08810          
       eta1              0.29048  0.11129    2.61001          
    --------------------------------------------------
    AIC= 1958.803 
    ||score||^2= 3.388644e-07

Measured by AIC there is a slight improvement in the model fit using the mixture model

    e0 <- estimate(m1,data=d)
    AIC(e0,em0)

        df      AIC
    e0  10 1961.839
    em0 12 1958.803

The spline model may then be estimated as before with the `two-stage` method

    em2 <- twostage(em0,m2,data=d)
    em2

                         Estimate Std. Error   Z-value   P-value
    Measurements:                                               
       y2~eta2            0.97823    0.16299   6.00179 1.952e-09
       y3~eta2            1.04530    0.10923   9.57015    <1e-12
    Regressions:                                                
       eta2~z             1.02885   20.08537   0.05122    0.9591
       eta2~eta1_1        2.80407   44.84864   0.06252    0.9501
       eta2~eta1_2       -0.02248    6.57588  -0.00342    0.9973
       eta2~eta1_3       -0.17334   19.50726  -0.00889    0.9929
       eta2~eta1_4        0.38674   22.12430   0.01748    0.9861
    Intercepts:                                                 
       y2                -0.12171    0.11148  -1.09173     0.275
       y3                -0.09870    0.10580  -0.93292    0.3509
       eta2               2.12361  112.75075   0.01883     0.985
    Residual Variances:                                         
       y1                 1.31872    0.77917   1.69247          
       y2                 1.09691    1.23383   0.88902          
       y3                 0.81345    1.13827   0.71464          
       eta2               1.99590    9.18558   0.21729

In this example the results are very similar to the Gaussian model:

    plot(I(eta2-z) ~ eta1, data=d, col=Col("black",0.5), pch=16,
         xlab=expression(eta[1]), ylab=expression(eta[2]))
    
    lines(Estimate ~ eta1, data=as.data.frame(p), col="darkblue", lwd=5)
    confband(p[,1], lower=p[,4], upper=p[,5], polygon=TRUE,
    	 border=NA, col=Col("darkblue",0.2))
    
    pm <- cbind(eta1=newd$eta1,
    	    estimate(em2, f=function(p) predict(e2,p=p,newdata=newd))$coefmat)
    lines(Estimate ~ eta1, data=as.data.frame(pm), col="darkred", lwd=5)
    confband(pm[,1], lower=pm[,4], upper=pm[,5], polygon=TRUE,
    	 border=NA, col=Col("darkred",0.2))
    legend("bottomright", c("Gaussian","Mixture"),
       col=c("darkblue","darkred"), lwd=2, bty="n")


<ref.bib>

