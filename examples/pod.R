library(stochbb)

#
# Processing model
#
d  <- 300
C  <- gamma(3,20)
X1 <- gamma(10,30)
# Shifted gamma
X2 <- gamma(1,70) %+% d


# response latency control condition
#   is simply R = C + X1
# The resulting distribution of response latecies is an analytic function on (0, infinity).
Rc <- C %+% X1

# And for the experimental condition 
#   R = C + min(X1, X2).
# This might be considered as a two parallel proceses X1 (identical to the X1 under control
# condition) and X2 which is delayed by d=300 ms but once started, much faster than X1.
# The resulting distribution is not analytic on (0, infinity)! It is identical to the distribution
# of R under control condition (Rc) on the interval (0, d] and diverges for T>d.
Re <- C %+% minimum(X1, X2)

Tmin <- 0; Tmax <- 1200; N <- 1200;
t   <- seq(Tmin, Tmax, length.out=N);
Tc  <- array(0, c(N)); Rc$density$eval(Tmin, Tmax, Tc)
Te  <- array(0, c(N)); Re$density$eval(Tmin, Tmax, Te)
TCc <- array(0, c(N)); Rc$density$evalCDF(Tmin, Tmax, TCc)
TCe <- array(0, c(N)); Re$density$evalCDF(Tmin, Tmax, TCe)

# Check with exact sampler
sam <- new(ExactSampler, c(Rc, Re))
res <- array(0, c(10000,2))
sam$sample(res)

# PDFs
plot(t, Te, type="l", lty=2)
lines(t, Tc, lty=1)
abline(v=d, lty=3)

# emp. PDFs (10000 sample)
plot(density(res[,2]), lty=2)
lines(density(res[,1]), lty=1)
abline(v=d, lty=3)

# Survival plots
plot(t, 1-TCc, type="l")
lines(t, 1-TCe, lty=2)
abline(v=d, lty=3)

# uncomment to plot "zoom" to divergence point.
#N <- 100000; dt <- (Tmax-Tmin)/N
#n0 <- (d-1-Tmin)/dt; n1 <- (d+1-Tmin)/dt
#t   <- seq(Tmin, Tmax, length.out=N);
#Tc  <- array(0, c(N)); Rc$density$eval(Tmin, Tmax, Tc)
#Te  <- array(0, c(N)); Re$density$eval(Tmin, Tmax, Te)
#plot(t[n0:n1], log(Te[n0:n1]/Tc[n0:n1]), type="l")
#abline(v=d, lty=3)
#abline(h=0, lty=4)
