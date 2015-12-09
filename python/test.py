from numpy import *
import stochbb

# Construct a normal(100, 30) RV
X = stochbb.normal(100, 30)
# eval its density
fX = empty(100,)
X.density().eval(fX, 0, 200)

# Create a compound-normal RV
Y = stochbb.normal(X, 30)
fY = empty(100,)
Y.density().eval(fY, 0, 200)

from matplotlib import pylab
t = linspace(0, 200, 200)
pylab.plot(t, fX)
pylab.plot(t, fY)
pylab.show()
