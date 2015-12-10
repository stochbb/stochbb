#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *
import stochbb
from matplotlib import pylab

#
# Global parameters
#

Tmin  = 0.;
Tmax  = 1500.;
steps = 1000;

# log word frequency (high, medium, low)
log_f = [3.0, 2., 1.]
# predictability (high, low, medium)
pred = [0.8, 0.01, 0.1]


#
# Early processing stages
#

# 90ms visual stage delay
# This however, I am not sure. Eq. 1 in Reichle 2003 makes absolutely no sense:
#  visual processing = t/(\epsilon ^ {\Simga i \_letter i - fixation\_/N })
# Oh boy! And this has been published.
visual = 90.


#
# lexical access
#
def lexical(start, log_f, pred):
    # first stage
    L1_mean = (228. - 10.*log_f)*(1-0.5*pred)
    L1_sd   = 0.18*L1_mean
    # Compute gamma distribution parameters
    L1_shape = (1/0.18)**2
    L1_scale = L1_mean/L1_shape
    # instantiate stage
    L1 = stochbb.gamma(L1_shape, L1_scale)

    # second stage
    L2_mean = 0.5*(228. - 10.*log_f)*(1-pred)
    L2_sd   = 0.18*L2_mean
    # Compute gamma distribution parameters
    L2_shape = (1/0.18)**2
    L2_scale = L2_mean/L2_shape
    # instantiate stage
    L2 = stochbb.gamma(L2_shape, L2_scale)
    # Retruns "events":
    #  first: "program saccade"
    #  second: "shift attention to next word"
    return (start+L1, start+L1+L2)


# lex processing of word N+0
X1, X12 = lexical(visual, log_f[0], pred[0])
# lex processing of word N+1 (triggered by attention shift)
X2, X22 = lexical(X12, log_f[1], pred[1])
# lex processing of word N+2 (triggered by attention shift)
X3, X32 = lexical(X22, log_f[1], pred[1])

#
# Motor control stages
#

# first stage: mean specified in paper, sd not
# what the heck? Anyway, assume the same sd like in L1 & L2
M1_mean = 187.
M1_sd   = 0.18*M1_mean
M1_shape = (1/0.18)**2;
M1_scale = M1_mean/M1_shape
# instantiate stages (3 i.i.d. realizations)
M11 = stochbb.gamma(M1_shape, M1_scale)
M12 = stochbb.gamma(M1_shape, M1_scale)
M13 = stochbb.gamma(M1_shape, M1_scale)

# motor control is a mixtrue in this case, here, there are at least
# 3 possible paths that can lead to a saccade. Hence we need to compute
# the probabilities for each case

# case 1; M11+X1 finsihed before X2, (M11+X1<X2) -> M11+X1+M2+...
# that is P((M11+X1)<X2) = \int F_{M11+X1}(t)\,f_{X3}(t)\,dt
C1 = M11+X1
F = empty(steps,); f = empty(steps,);
C1.density().evalCDF(Tmin, Tmax, F);
X2.density().eval(Tmin, Tmax, f);
p1 = sum(F*f);
# case 2: (M11+X1>X2) & (M12+X2<X3) -> M12+X2+M2+...
C2 = M12+X2
C2.density().evalCDF(Tmin, Tmax, F);
X3.density().eval(Tmin, Tmax, f);
p2 = (1-p1)*sum(F*f)
# case 3: (M11+X1>X2) & (M12+X2>X3) -> M13+X3+M2+...
C3 = M13+X3
p3 = (1-p1)*(1-sum(F*f))
print( (p1, p2, p3) )

# second stage
M2_mean  = 53.
M2_sd    = 0.18*M2_mean
M2_shape = (1/0.18)**2;
M2_scale = M2_mean/M2_shape
M2 = stochbb.gamma(M2_shape, M2_scale)
# and finally saccade gen. (fixed 25ms)
S = 25.

# Now assemble fixation duration RV as the mixture
dur = S + M2 + stochbb.mixture(p1,C1, p2,C2, p3,C3)
# Only the case w/o attention shift
path1 = S + M2 + C1
# Only the case w/ attention shift to N+1
path2 = S + M2 + C2
# Only the case w/ attention shift to N+1 and then to N+2
path3 = S + M2 + C3

t = linspace(Tmin, Tmax, steps)
# Eval density of mixture
pdf = empty(steps,);    dur.density().eval(Tmin, Tmax, pdf)
# eval density w/o attention shift
p1_pdf = empty(steps,); path1.density().eval(Tmin, Tmax, p1_pdf);
# eval density w/ attention shift to N+1
p2_pdf = empty(steps,); path2.density().eval(Tmin, Tmax, p2_pdf);
# eval density w/ attention shift to N+1 and then to N+2
p3_pdf = empty(steps,); path3.density().eval(Tmin, Tmax, p3_pdf);

pylab.plot(t, pdf)
pylab.plot(t, p1_pdf)
pylab.plot(t, p2_pdf)
pylab.plot(t, p3_pdf)
pylab.show()
