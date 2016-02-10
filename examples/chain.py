#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import stochbb;
from matplotlib import pylab

#
# Corpus model
#

# Log word frequency (per mililon)
#  is uniform distr. on the interval [0,4]
f = stochbb.uniform(0,4)

# Predictability is unif. distributed on the interval [0,1]
p = stochbb.uniform(0,1)

#
# Cognitive model
#

# "lexical" stage is gamma distributed with k=3*f+5, theta=10
L = stochbb.gamma(3*f+5, 10);
# "semantic" stage is gamma distributed with k=3*p+5, theta=20
S = stochbb.gamma(5*p+5, 20);
# motor stage, gamma distributed with k=10, theta=10
M = stochbb.gamma(10,10)

# response latency is simply R = L + S + M
R = L+S+M

Tmin, Tmax, N = 0, 500, 10000;
t = linspace(Tmin, Tmax, N); dt = float(Tmax-Tmin)/N
pL = empty(N,); L.density().eval(Tmin, Tmax, pL)
pS = empty(N,); S.density().eval(Tmin, Tmax, pS)
pM = empty(N,); M.density().eval(Tmin, Tmax, pM)
pR = empty(N,); R.density().eval(Tmin, Tmax, pR)


pylab.plot(t, pL, label="Lexical")
pylab.plot(t, pS, label="Semantic")
pylab.plot(t, pM, label="Motor")
pylab.plot(t, pR, label="Response (L+S+M)")
pylab.legend()
pylab.show()
