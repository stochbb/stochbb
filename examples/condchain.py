#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib import pylab
import stochbb

X1 = stochbb.gamma(3, 100);
X2 = stochbb.gamma(3, 120);
Y1 = stochbb.gamma(3, 100);

Z1 = stochbb.condchain(X1, X2, Y1, Y1)
Z2 = stochbb.minimum(X1,X2)+Y1

Tmin, Tmax, steps = 0., 1200., 1000
t = linspace(Tmin, Tmax, steps)
p1 = empty(steps,); Z1.density().eval(Tmin, Tmax, p1);
p2 = empty(steps,); Z2.density().eval(Tmin, Tmax, p2);

pylab.plot(t, p1)
pylab.plot(t, p2)
pylab.show()
