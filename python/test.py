from numpy import *
import stochbb

# Construct a normal(100, 30) RV
t = linspace(0, 200, 1000)
X = stochbb.normal(50, 10)
Y = stochbb.normal(X, 10)
fY = empty(1000,)
Y.density().eval(0, 200, fY)

from matplotlib import pylab
pylab.plot(t, fZ)
pylab.show()
