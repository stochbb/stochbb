library(stochbb)
Tmin <- 0; Tmax <- 700; N <- 1000;
t = seq(Tmin, Tmax, length.out=N)

X <- gamma(10,10);
Y <- gamma(10,20);
Z = X %+% Y %+% 100;

rX <- X$density$rangeEst(0.01);
rY <- Y$density$rangeEst(0.01);
rZ <- Z$density$rangeEst(0.01);

dX <- array(0, c(N)); X$density$eval(Tmin, Tmax, dX);
dY <- array(0, c(N)); Y$density$eval(Tmin, Tmax, dY);
dZ <- array(0, c(N)); Z$density$eval(Tmin, Tmax, dZ);

plot(t, dX, type="l"); abline(v=rX, lty=1)
lines(t, dY, lty=2); abline(v=rY, lty=2)
lines(t, dZ, lty=3); abline(v=rZ, lty=3)

