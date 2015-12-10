#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *
import stochbb
from matplotlib import pylab

#
# Global parameters
#

# log word frequency
log_f = [2.0, 2.3, 1.9]
# predictability
pred = [0.1, 0.2, 0.06]


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

#
# Motor control
#

# first stage, mean specified in paper, sd not
# what the heck? Anyway, assume the same shit like in L1 & L2
M1_mean = 187.
M1_sd   = 0.18*M1_mean
M1_shape = (1/0.18)**2;
M1_scale = M1_mean/M1_shape
# instantiate stages (3 indp realizations)
M11 = stochbb.gamma(M1_shape, M1_scale)
M12 = stochbb.gamma(M1_shape, M1_scale)
M13 = stochbb.gamma(M1_shape, M1_scale)

# second stage
M2_mean  = 53.
M2_sd    = 0.18*M2_mean
M2_shape = (1/0.18)**2;
M2_scale = M2_mean/M2_shape
M2 = stochbb.gamma(M2_shape, M2_scale)

# saccade gen. (fixed 25ms)
S = 25.


# lex processing of word 0
X1, X12 = lexical(visual, log_f[0], pred[0])
# lex processing of word 1 (trigger by attention shift)
X2, X22 = lexical(X12, log_f[1], pred[1])
# lex processing of word 2 (trigger by attention shift)
X3, X32 = lexical(X22, log_f[1], pred[1])

# motor control (what ever comes first)
dur = S + M2 + M13 + stochbb.minimum(M12 + stochbb.minimum(M11 + X1, X2), X3)
