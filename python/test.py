from numpy import *
import stochbb

# Construct a normal(100, 30) RV
X = stochbb.normal(100, 30)
sampler = stochbb.MarginalSampler(X, 0, 200, 1000)
samples = empty(1000,)
sampler.sample(samples)
print(samples)

fX = empty(1000,)
X.density().eval(0, 200, fX)


