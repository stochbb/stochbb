#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib import pylab
import stochbb

X1 = stochbb.gamma(3, 100) + 100;
X2 = stochbb.gamma(3, 120) + 100;
Y1 = stochbb.gamma(3, 140) + 100;
Y2 = stochbb.gamma(3, 130) + 100;

Z1 = stochbb.condchain(X1, X2, Y1, Y2)

Tmin, Tmax, steps = 0, 2000., 100
t = linspace(Tmin, Tmax, steps)

pX1 = zeros(steps,); X1.density().eval(Tmin, Tmax, pX1);
pX2 = zeros(steps,); X2.density().eval(Tmin, Tmax, pX2);
pY1 = zeros(steps,); Y1.density().eval(Tmin, Tmax, pY1);
pY2 = zeros(steps,); Y2.density().eval(Tmin, Tmax, pY2);
pZ1 = zeros(steps,); Z1.density().eval(Tmin, Tmax, pZ1);

pylab.plot(t, pX1, label="X1")
pylab.plot(t, pX2, label="X2")
pylab.plot(t, pY1, label="Y1")
pylab.plot(t, pY2, label="Y2")
pylab.plot(t, pZ1, label="Z1")
pylab.legend()

pylab.show()
