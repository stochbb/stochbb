#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib import pylab
import stochbb

X1 = stochbb.gamma(3, 100);
X2 = stochbb.gamma(3, 120);
Y1 = stochbb.gamma(3, 140);

# this is equivalent to Y1 + min(X1,X2)
Z1 = stochbb.condchain(X1, X2, Y1, Y1)
# compare
Z2 = stochbb.minimum(X1,X2)+Y1
Z3 = stochbb.minimum(X1,X2)

print Z1.density()
print Z2.density()

Tmin, Tmax, steps = 0, 1500., 100
t = linspace(Tmin, Tmax, steps)
pX1 = zeros(steps,); X1.density().eval(Tmin, Tmax, pX1);
pX2 = zeros(steps,); X2.density().eval(Tmin, Tmax, pX2);
pY1 = zeros(steps,); Y1.density().eval(Tmin, Tmax, pY1);
pZ1 = zeros(steps,); Z1.density().eval(Tmin, Tmax, pZ1);
pZ2 = zeros(steps,); Z2.density().eval(Tmin, Tmax, pZ2);

pylab.subplot(211)
pylab.plot(t, pX1, label="X1")
pylab.plot(t, pX2, label="X2")
pylab.plot(t, pY1, label="Y1")
pylab.plot(t, pZ1, label="Z1")
pylab.plot(t, pZ2, label="Z2")
pylab.legend()


Tmin, Tmax, steps = -200, 1500., 100
t = linspace(Tmin, Tmax, steps)
pX1 = zeros(steps,); X1.density().eval(Tmin, Tmax, pX1);
pX2 = zeros(steps,); X2.density().eval(Tmin, Tmax, pX2);
pY1 = zeros(steps,); Y1.density().eval(Tmin, Tmax, pY1);
pZ1 = zeros(steps,); Z1.density().eval(Tmin, Tmax, pZ1);
pZ2 = zeros(steps,); Z2.density().eval(Tmin, Tmax, pZ2);

pylab.subplot(212)
pylab.plot(t, pX1, label="X1")
pylab.plot(t, pX2, label="X2")
pylab.plot(t, pY1, label="Y1")
pylab.plot(t, pZ1, label="Z1")
pylab.plot(t, pZ2, label="Z2")
pylab.legend()

pylab.show()
