from numpy import *
import stochbb
from matplotlib import pylab

X1 = stochbb.gamma(2,1.1)+1
X2 = stochbb.gamma(2,1)+1
X3 = X1 + X2

Tmin, Tmax, N = 0, 4, 100
t = linspace(Tmin, Tmax, N)
dX1 = empty((N,)); X1.density().eval(Tmin, Tmax, dX1)
dX2 = empty((N,)); X2.density().eval(Tmin, Tmax, dX2)
dX3 = empty((N,)); X3.density().eval(Tmin, Tmax, dX3)

pylab.plot(t, dX1, label="X1")
pylab.plot(t, dX2, label="X2")
pylab.plot(t, dX3, label="X1+X2")


Tmin, Tmax, N = 0, 12, 300
t = linspace(Tmin, Tmax, N)
dX12 = empty((N,)); X1.density().eval(Tmin, Tmax, dX12)
dX22 = empty((N,)); X2.density().eval(Tmin, Tmax, dX22)
dX32 = empty((N,)); X3.density().eval(Tmin, Tmax, dX32)

pylab.plot(t, dX12, "--", label="X1")
pylab.plot(t, dX22, "--", label="X2")
pylab.plot(t, dX32, "--", label="X1+X2")

pylab.legend()
pylab.show()

