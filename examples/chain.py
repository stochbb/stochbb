#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
import stochbb;
from matplotlib import pylab

X1 = stochbb.gamma(10, 10);
X2 = stochbb.gamma(10, 20);
Y = X1+X2

Tmin, Tmax, N = 0, 500, 10000;
t = linspace(Tmin, Tmax, N); dt = float(Tmax-Tmin)/N
pX1 = empty(N,); X1.density().eval(Tmin, Tmax, pX1)
pX2 = empty(N,); X2.density().eval(Tmin, Tmax, pX2)
pY = empty(N,); Y.density().eval(Tmin, Tmax, pY)


pylab.plot(t, pX1, label="X1~Gamma(10,10)")
pylab.plot(t, pX2, label="X2~Gamma(10,20)")
pylab.plot(t, pY, label="Y=X1+X2")
pylab.legend()
pylab.show()
