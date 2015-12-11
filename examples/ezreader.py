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
dt = float(Tmax-Tmin)/steps

# log word frequency (high, medium, low)
log_f = [3.0, 2., 1.]
# predictability (high, high, medium)
pred = [0.8, 0.6, 0.1]


#
# Early processing stages
#

# 90ms visual stage delay
# This however, I am not sure. Eq. 1 in Reichle 2003 makes absolutely no sense:
#  visual processing = t/(\epsilon ^ {\Simga i \_letter i - fixation\_/N })
# Oh boy! And this has been published, where t is the fixation duration.
visual = 90.


#
# lexical access variables
#  this function constructs the L1 and L2 RVs given the frequency and predictability
#  of the attended word
def lexical(log_f, pred):
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
    return (L1, L2)


# lex processing of word N+0
L11, L21 = lexical(log_f[0], pred[0])
# lex processing of word N+1 (triggered by attention shift)
L12, L22 = lexical(log_f[1], pred[1])
# lex processing of word N+2 (triggered by attention shift)
L13, L23 = lexical(log_f[2], pred[2])

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

#
# The motor control is a mixtrue in this case. Here, there are at least
# 3 possible paths that can lead to a saccade.
# ... -> L11 -> M11 -> M2 -> ...,                              : No attention shift
# ... -> L11 -> L21 -> L12 -> M12 -> M2 -> ... and             : attention shift to N+1
# ... -> L11 -> L21 -> L12 -> L22 -> L13 -> M13 -> M2 -> ...   : attention shift to N+2
#
# case 1; M11 finsihed before L21 & L12, (M11<L21+L12)
# case 2: (M11>L12+L21) & (M12<L22+L13)
# case 3: (M11>L12+L21) & (M12>L22+L13)
P1 = visual + L11 + M11;
P2 = visual + L11 + L21 + L12 + M12;
P3 = visual + L11 + L21 + L12 + L22 + L13 + M13


# second stage
M2_mean  = 53.
M2_sd    = 0.18*M2_mean
M2_shape = (1/0.18)**2;
M2_scale = M2_mean/M2_shape
M2 = stochbb.gamma(M2_shape, M2_scale)
# and finally saccade gen. (fixed 25ms)
S = 25.

# Now assemble fixation duration RV as the mixture
dur = S + M2 + stochbb.conditional(M11, L21+L12, P1, stochbb.conditional(M12, L22+L13, P2, P3));

# update paths
P1 += S + M2;
P2 += S + M2;
P3 += S + M2;

t = linspace(Tmin, Tmax, steps)
# Eval density of mixture
pdf = empty(steps,);    dur.density().eval(Tmin, Tmax, pdf)
# eval density w/o attention shift
P1_pdf = empty(steps,); P1.density().eval(Tmin, Tmax, P1_pdf);
# eval density w/ attention shift to N+1
P2_pdf = empty(steps,); P2.density().eval(Tmin, Tmax, P2_pdf);
# eval density w/ attention shift to N+1 and then to N+2
P3_pdf = empty(steps,); P3.density().eval(Tmin, Tmax, P3_pdf);

pylab.plot(t, pdf, label="EZ-Reader")
pylab.plot(t, P1_pdf, label="no shift")
pylab.plot(t, P2_pdf, label="shift to N+1")
pylab.plot(t, P3_pdf, label="shift to N+2")
pylab.xlabel("Time [ms]")
pylab.ylabel("f(T)")
pylab.legend()
pylab.savefig("ezreader.pdf")
pylab.show()
