library(stochbb)

# Var + x (x numeric) is eqiv. to affine(Var, 1, x)
X1 <- gamma(3, 100) %+% 100;
X2 <- gamma(3, 120) %+% 100;
Y1 <- gamma(3, 140) %+% 100;
Y2 <- gamma(3, 130) %+% 100;

# Z1 = Y1+X1 if X1 < X2, Y2+X2 else.
Z1 <- condchain(X1, X2, Y1, Y2)

Tmin <-0; Tmax<-2000; steps <- 100;
t = seq(Tmin, Tmax, length.out=steps)

pX1 <- array(0, c(steps)); X1$density$eval(Tmin, Tmax, pX1);
pX2 <- array(0, c(steps)); X2$density$eval(Tmin, Tmax, pX2);
pY1 <- array(0, c(steps)); Y1$density$eval(Tmin, Tmax, pY1);
pY2 <- array(0, c(steps)); Y2$density$eval(Tmin, Tmax, pY2);
pZ1 <- array(0, c(steps)); Z1$density$eval(Tmin, Tmax, pZ1);

plot(t, pX1, type="l")
lines(t, pX2, lty=1)
lines(t, pY1, lty=2)
lines(t, pY2, lty=3)
lines(t, pZ1, lty=4)
