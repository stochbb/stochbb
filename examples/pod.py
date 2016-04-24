#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import stochbb;
from matplotlib import pylab
import cProfile

# Install Logger
stochbb.Logger.addHandler(stochbb.IOLogHandler())

#
# Processing model
#
d = 300
C = stochbb.gamma(5,20)
X1 = stochbb.gamma(10,30)
D  = stochbb.delta(d)
X2 = stochbb.gamma(3,50)


# response latency control condition
#   is simply R = C + X1
Rc = C+X1
# and for the experimental condition
Re = C+stochbb.minimum(X1, D+X2)

Tmin, Tmax, N = 0, 1200, 1200;
t = linspace(Tmin, Tmax, N);
Tc = empty(N,); Rc.density().eval(Tmin, Tmax, Tc)
Te = empty(N,); Re.density().eval(Tmin, Tmax, Te)
TCc = empty(N,); Rc.density().evalCDF(Tmin, Tmax, TCc)
TCe = empty(N,); Re.density().evalCDF(Tmin, Tmax, TCe)


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
