#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import stochbb;
from matplotlib import pylab
import cProfile

# Install Logger
#stochbb.Logger.addHandler(stochbb.IOLogHandler())

#
# Corpus model
#

# Log word frequency (per mililon)
#  is uniform distr. on the interval [0,4]
f = stochbb.uniform(0,4)

# Predictability is unif. distributed on the interval [0,1]
#  somewhat unrealistic but this is just an example
p = stochbb.uniform(0,1)

#
# Cognitive model
#

# "lexical" stage is gamma distributed with k=5*f+5, theta=10
L = stochbb.gamma(5*f+5, 10);
# "semantic" stage is gamma distributed with k=10*p+5, theta=20
S = stochbb.gamma(10*p+5, 20);
# motor stage, gamma distributed with k=10, theta=30
M = stochbb.gamma(10, 30)

# response latency is simply R = L + S + M
R = L+S+M

print R, R.density()
Tmin, Tmax, N = 0, 1200, 1000;
t = linspace(Tmin, Tmax, N); dt = float(Tmax-Tmin)/N
pL = empty(N,); L.density().eval(Tmin, Tmax, pL)
pS = empty(N,); S.density().eval(Tmin, Tmax, pS)
pM = empty(N,); M.density().eval(Tmin, Tmax, pM)
pR = empty(N,); cProfile.run("R.density().eval(Tmin, Tmax, pR)")


X = empty((100000, 4))
sam = stochbb.ExactSampler([L,S,M,R])
cProfile.run("sam.sample(X)")

ll, = pylab.plot(t, pL, "b", lw=3)
pylab.hist(X[:,0], color="b", bins=100, normed=True)
ls, = pylab.plot(t, pS, "g", lw=3)
pylab.hist(X[:,1], color="g", bins=100, normed=True)
lm, = pylab.plot(t, pM, "r", lw=3)
pylab.hist(X[:,2], color="r", bins=100, normed=True)
lr, = pylab.plot(t, pR, "c", lw=3)
pylab.hist(X[:,3], color="c", bins=100, normed=True)
pylab.legend((ll, ls, lm, lr),
             ("Lexical", "Semantic", "Motor", "Response (L+S+M)"))

pylab.savefig("example.pdf")
pylab.show()
