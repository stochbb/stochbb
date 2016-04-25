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

# "lexical" stage is gamma distributed with k=5*f+5, theta=10
L <- gamma(affine(f, 5, 5), 10);
# "semantic" stage is gamma distributed with k=10*p+5, theta=20
S <- gamma(affine(p, 10, 5), 20);
# motor stage, gamma distributed with k=10, theta=30
M <- gamma(10, 30)

# response latency is simply R = L + S + M
R = L %+% S %+% M

Tmin <- 0; Tmax <- 1200; N <- 1000;
t  <- seq(Tmin, Tmax, length.out=N);
dt <- (Tmax-Tmin)/N

pL <- rep(0,N);  L$density$eval(Tmin, Tmax, pL)
pS <- rep(0,N);  S$density$eval(Tmin, Tmax, pS)
pM <- rep(0,N);  M$density$eval(Tmin, Tmax, pM)
pR <- rep(0,Ns); R$density$eval(Tmin, Tmax, pR)


#X = empty((100000, 4))
#sam = stochbb.ExactSampler([L,S,M,R])

plot(t, pL, type="L")
#pylab.hist(X[:,0], color="b", bins=100, normed=True)
lines(t, pS, lty=2)
#pylab.hist(X[:,1], color="g", bins=100, normed=True)
lines(t, pM, lty=3)
#pylab.hist(X[:,2], color="r", bins=100, normed=True)
lines(t, pR, lty=4)
#pylab.hist(X[:,3], color="c", bins=100, normed=True)
