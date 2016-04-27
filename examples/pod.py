#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import stochbb;
from matplotlib import pylab
import cProfile

#
# This example demonstrates how a divergence point can arise in a very simple stochastic
# model of some cognitive process. The divergence point is the earliest time-point at which the
# distributions of some response latency differ between two experimental conditions. Obviously,
# there is no divergence point if the two response latency distributions are analytic. For two
# non-analytic response latencies, a divergence point may exist.
#
# This example uses 3 Gamma distributions (of wich one is a shifted Gamma distribution) and the
# minimum(,) operator. For both conditions (control, experimental) there is a common stage C
# following some Gamma distribution. Under control condition, the common stage is followed by a
# second Gamma distributed stage X1 (C -> X1). Under the experimental condition, the common stage
# also triggers a second stage X2 with a shifted Gamma distribution parallel with the stage X1
# present under the control conditon. This stage is delayed by 300ms relative to the parallel stage
# X1, but completes much faster than X1 once it started. The response latency under the exp.
# condition is then the chain C -> max(X1, X2).
#

#
# Processing model
#
d = 300
C = stochbb.gamma(3,20)
X1 = stochbb.gamma(10,30)
# shifted gamma
X2 = stochbb.gamma(1,70) + d


# response latency control condition
#   is simply R = C + X1
Rc = C + X1
# and for the experimental condition
Re = C + stochbb.minimum(X1, X2)

Tmin, Tmax, N = 0, 1200, 1200;
t = linspace(Tmin, Tmax, N);
Tc = empty((N,)); Rc.density().eval(Tmin, Tmax, Tc)
Te = empty((N,)); Re.density().eval(Tmin, Tmax, Te)
TCc = empty((N,)); Rc.density().evalCDF(Tmin, Tmax, TCc)
TCe = empty((N,)); Re.density().evalCDF(Tmin, Tmax, TCe)

pylab.subplot(211)
pc, = pylab.plot(t, Tc, "b", lw=3)
pe, = pylab.plot(t, Te, "g", lw=3)
pod = pylab.axvline(d, color="r")
pylab.legend((pc, pe, pod),
             ("Control", "Experimental", "PoD"))

pylab.subplot(212)
pc, = pylab.plot(t, 1-TCc, "b", lw=3)
pe, = pylab.plot(t, 1-TCe, "g", lw=3)
pod = pylab.axvline(d, color="r")
pylab.legend((pc, pe, pod),
             ("Control", "Experimental", "PoD"))
pylab.savefig("pod.pdf")
pylab.show()
