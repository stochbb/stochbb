from numpy import *
from matplotlib import pylab
import stochbb

X1 = stochbb.normal(100,10)
X2 = stochbb.normal(100,10)
# Determine res numerically
d  = stochbb.directConvolve(X1.density(), X2.density())
# Will determine density analytically
Y = X1+X2
# What should be the solution
Z = stochbb.normal(200, sqrt(200))

Tmin, Tmax, N = 0, 400, 4000
t = linspace(Tmin, Tmax, N)

f1 = empty((N,)); d.eval(Tmin, Tmax, f1)
f2 = empty((N,)); Y.density().eval(Tmin, Tmax, f2)
f3 = empty((N,)); Z.density().eval(Tmin, Tmax, f3)

p1, = pylab.plot(t, f1)
p2, = pylab.plot(t, f2)
p3, = pylab.plot(t, f3)
pylab.legend((p1, p2, p3),
             ("Numerical", "Analytic (internal)", "Analytic (true)"))
pylab.show()
