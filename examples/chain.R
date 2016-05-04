library(stochbb)

#
# Corpus model
#

# Log word frequency (per mililon)
#  is uniform distr. on the interval [0,4]
f <- uniform(0,4)

# Predictability is unif. distributed on the interval [0,1]
#  somewhat unrealistic but this is just an example
p <- uniform(0,1)

#
# Cognitive model
#

# "lexical" stage is compound gamma distributed with k=5*f+5, theta=10
L <- gamma(affine(f, 5, 5), 10);
# "semantic" stage is compound gamma distributed with k=10*p+5, theta=20
S <- gamma(affine(p, 10, 5), 20);
# motor stage, simple gamma distributed with k=10, theta=30
M <- gamma(10, 30)

# response latency is simply R = L + S + M
# can also be expressed as R <- L %+% S %+% M
R <- chain(L, S, M)

Tmin <- 0; Tmax <- 1200; N <- 1000;
t  <- seq(Tmin, Tmax, length.out=N);

pL <- array(0, c(N)); L$density$eval(Tmin, Tmax, pL)
pS <- array(0, c(N)); S$density$eval(Tmin, Tmax, pS)
pM <- array(0, c(N)); M$density$eval(Tmin, Tmax, pM)
pR <- array(0, c(N)); R$density$eval(Tmin, Tmax, pR)


sam <- new(ExactSampler, list(L, S, M, R))
X <- matrix(0, 100000, 4); sam$sample(X)

plot(t, pL, type="l")
lines(density(X[,1]), lty=1)
lines(t, pS, lty=2)
lines(density(X[,2]), lty=2)
lines(t, pM, lty=3)
lines(density(X[,3]), lty=3)
lines(t, pR, lty=4)
lines(density(X[,4]), lty=4)
