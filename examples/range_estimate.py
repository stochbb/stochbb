from numpy import *
import stochbb
from matplotlib import pylab

Tmin, Tmax, N = 0, 700, 1000;
t = linspace(Tmin, Tmax, N)

stochbb.Logger.addHandler(stochbb.IOLogHandler())

X = stochbb.gamma(10,10);
Y = stochbb.gamma(10,20);
Z = X+Y+100;
print X.density().rangeEst(0.01)
print Y.density().rangeEst(0.01)
print Z.density().rangeEst(0.01)

dX = empty(N,); X.density().eval(Tmin, Tmax, dX);
dY = empty(N,); Y.density().eval(Tmin, Tmax, dY);
dZ = empty(N,); Z.density().eval(Tmin, Tmax, dZ);

pylab.plot(t, dX)
pylab.plot(t, dY)
pylab.plot(t, dZ)

pylab.show()

