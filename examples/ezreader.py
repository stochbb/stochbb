#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import *
import stochbb
import matplotlib
matplotlib.use("Agg")
from matplotlib import pylab

#
# Global parameters
#

Tmin  = 0.0;
Tmax  = 1500.;
steps = 1000;


def ezreader(log_f, pred):
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
    #  of the words
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

        # done
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

    # second stage
    M2_mean  = 53.
    M2_sd    = 0.18*M2_mean
    M2_shape = (1/0.18)**2;
    M2_scale = M2_mean/M2_shape
    M2 = stochbb.gamma(M2_shape, M2_scale)
    # and finally saccade gen. (fixed 25ms)
    S = 25.

    #
    # The motor control is a mixtrue in this case. Here, there are at least
    # 3 possible paths that can lead to a saccade.
    #
    # case 1; M11 finsihed before L21 & L12, (M11<L21+L12)
    # ... -> M11 -> ...,                              : No attention shift
    P1 = visual + L11 + M11 + (M2 + S);
    # case 2: (M11>L12+L21) & (M12<L22+L13)
    # ... -> L21 -> L12 -> M12 -> ... and             : attention shift to N+1
    P2 = visual + L11 + L21 + L12 + M12 + (M2 + S);
    # case 3: (M11>L12+L21) & (M12>L22+L13)
    # ... -> L21 -> L12 -> L22 -> L13 -> M13 -> ...   : attention shift to N+2
    P3 = visual + L11 + L21 + L12 + L22 + L13 + M13 + (M2 + S);
    # case 4: You could also add a 3rd attention shift here but I think that it is
    #         unlikely

    # Now assemble fixation duration RV as a conditional chained RV
    dur = visual + L11;
    dur += stochbb.condchain(M11, L21+L12,
                             stochbb.delta(0),
                             stochbb.condchain(M12, L22+L13,
                                               stochbb.delta(0),
                                               M13));
    dur += M2 + S
    # Done
    return dur;

# log word frequency (high, medium, low)
log_f = [3.0, 2., 1.]
# predictability (high, high, medium)
pred = [0.7, 0.6, 0.2]
# scan p2 (pred[1]) from 0 to 0.9 in 100 steps
p2_range = linspace(0, 0.9, 100)
# time scales
t = linspace(Tmin, Tmax, steps)

def eval_pdf(p):
    pred[1] = p;
    dur = ezreader(log_f, pred)
    # Eval density of mixture
    pdf = empty(steps,); dur.density().eval(Tmin, Tmax, pdf)
    return pdf;

fig, ax = pylab.subplots()
line, = ax.plot(t, eval_pdf(0))
text = ax.text(1000, 0.008, "p(N+1)=0")

def animate(p):
    line.set_ydata(eval_pdf(p));
    text.set_text("p(N+1)={0:.2f}".format(p))
    return line,text

def init():
    return animate(0);

import matplotlib.animation as animation
ani = animation.FuncAnimation(fig, animate, p2_range, init_func=init,
                              interval=25);

# Set up formatting for the movie files
Writer = animation.writers['avconv']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

ani.save("ezreader.mp4", writer=writer)
#pylab.show()
