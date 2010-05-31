###
### All lava-demos:
###

m <- lvm(y~x)
d <- sim(m,100)

e <- estimate(y~x,d)
e

