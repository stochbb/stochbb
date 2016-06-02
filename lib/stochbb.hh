/** @mainpage
 * StochBB provides a set of classes and functions (see @ref api) that allows
 * to assemble complex systems of random variables (RVs) and to obtain their PDFs and CDFs
 * numerically rather than by means of stochastic simulation. Additionally, it is possible to
 * sample from these RVs for verification. For some in-depth description of StochBB, consider
 * reading the <a href="https://hmatuschek.github.io/stochbb/manual/manual.pdf">manual</a>.
 *
 * Frequently, complex systems are described in terms of stochastic processes, as the
 * underlying deterministic process is too complex to be modeled exactly or as the
 * process is indeed random. It is not always the random process itself that is of
 * interest, but a derived quantity. For example, the distribution of waiting times until the
 * process reaches a certain state. In the field of cognitive psychology, random processes are
 * frequently used to describe each processing stage in a chain of stages that leads to a response.
 * The state of each random process itself is usually not measurable but the total response time of
 * all processing stages involved. Although each processing stage is modeled as a random process,
 * the waiting-time of a single stage is just a random variable and the complete system is then a
 * system of dependent random variables.
 *
 * StochBB is able to describe and analyze complex systems of dependent random variables
 * by combining simple ones (representing single stages with a known waiting-time distribution) to
 * a complex system. For example, consider the following independent random variables
 * \f[
 *   X_1 \sim \Gamma(10, 100)\,, X_2 \sim \Gamma(20, 50)\text{ and } X_3 \sim \text{Exp}(0.01)\,,
 * \f]
 * which are described completely by their distribution. This means that the time, the
 * processing stage \f$X_1\f$ needs to complete is gamma-distributed with shape \f$k=10\f$ and scale
 * \f$\theta=100\f$. Analogously, the stages \f$X_2\f$ and \f$X_3\f$ are defined by their own
 * waiting-time distribution.
 *
 * From these basic building blocks, a more complex system can be assembled by
 * combining them. Using the example above, one may define a new processing stage that is a chain of
 * the stages \f$X_1,\dots,X_3\f$. This simple chain then describes the successive
 * processing of information entering the first stage described by \f$X_1\f$.
 * Once the first processing stage finished, its result gets forwarded to the second stage represented
 * by \f$X_2\f$ and finally to the last stage represented by \f$X_3\f$. The waiting time of the
 * complete chain is again a random variable that is the sum of all random variables,
 * or expressed mathematically
 * \f[
 *   Y = X_1 + X_2 + X_3\,.
 * \f]
 *
 * SochBB determines the probability density function (PDF) or cumulative probability function
 * (CDF) of the waiting-time distribution of the process \f$Y\f$ analytically (as far as
 * possible) or resorts to a numeric method if the analytic approach fails. More over it provides
 * an efficient and correct sampler for the system of random variables.
 *
 * Continuing the example above, please note that the sum of random variables commutes. Hence the
 * random variable \f$Y\f$ remains the same if defined as \f$Y = X_1 + X_3 + X_2\f$ instead of
 * \f$Y = X_1 + X_2 + X_3\f$. Moreover, the distribution of the sum of \f$X_1\f$ and \f$X_3\f$ can
 * be determined analytically as \f$X' = (X_1+X_3)\sim\Gamma(11,100)\f$. Hence the random variable
 * \f$Y\f$ can now be expressed as \f$Y = X_2 + X'\f$, and only a single numeric convolution is
 * necessary to obtain the PDF of the random variable \f$Y\f$.
 *
 * StochBB implements several reductions of the system of random variables,
 * exploiting mathematical identities of random variables. To this end, it allows to obtain the
 * PDFs and CDFs of random variables efficiently.
 *
 * @author Hannes Matuschek <hannes [dot] matuschek [at] uni-potsdam [dot] de>
 */

#ifndef __SBB_STOCHBB_HH__
#define __SBB_STOCHBB_HH__

#include "logger.hh"
#include "api.hh"

#include "density.hh"
#include "randomvariable.hh"
#include "affinetrafo.hh"
#include "chain.hh"
#include "minmax.hh"
#include "mixture.hh"
#include "conditional.hh"
#include "compound.hh"
#include "exactsampler.hh"
#include "marginalsampler.hh"

#include "math.hh"

#endif // __SBB_STOCHBB_HH__

