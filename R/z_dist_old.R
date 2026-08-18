##' @export
coxWeibull.lvm <- dist_cox_weibull

##' @export
coxExponential.lvm <- dist_cox_exponential

##' @export
aalenExponential.lvm <- dist_aalen_exponential

##' @export
coxGompertz.lvm <- dist_cox_gompertz

##' @export
normal.lvm <- dist_gaussian

##' @export
gaussian.lvm <- normal.lvm

##' @export
lognormal.lvm <- dist_lognormal

##' @export
poisson.lvm <- dist_poisson

##' @export
pareto.lvm <- dist_pareto

##' @export
threshold.lvm <- dist_threshold

##' @export
multinomial.lvm <- dist_multinomial

##' @export
binomial.lvm <- dist_bernoulli

##' @export
logit.lvm <- binomial.lvm("logit")

##' @export
probit.lvm <- binomial.lvm("probit")

##' @export
Gamma.lvm <- dist_gamma

##' @export
loggamma.lvm <- dist_loggamma

##' @export
chisq.lvm <- dist_chisq

##' @export
student.lvm <- dist_t

##' @export
uniform.lvm <- dist_uniform

##' @export
weibull.lvm <- dist_weibull

##' @export
id.lvm <- dist_seqint

##' @export
Sequence.lvm <- dist_seq

##' @export
none.lvm <- dist_none

##' @export
constant.lvm <- dist_const

##' @export
ones.lvm <- dist_ones

##' @export
binary.lvm <- ones.lvm

##' @export
beta.lvm <- dist_beta

##' @export
mvn.lvm <- dist_mvn

##' @export
GM2.lvm <- dist_gaussian_mixture2

##' @export
GM3.lvm <- dist_gaussian_mixture3
