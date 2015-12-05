/** @mainpage Introduction
 * Frequently, complex systems are described in terms of stochastic processes. Either as the
 * underlaying deterministic process is simply too complex to be modelled exactly or as the
 * process is indeed random. The most simple class of random processes are the stationay processes.
 * There the process parameters do not change in time and the process can be described as a simple
 * random variable. StochBB can be used to describe and analyze complex stationary random processes
 * in time by combining simple processes (described by a known waiting-time distribution) to
 * complex processes.
 *
 * For example, consider the following indepdenent and stationary random processes
 * \f[
 *  X_1 \sim \Gamma(10, 100)\,,\quad X_2 \sim \Gamma(20, 50)\,, X_3 \sim \text{Exp}(0.01);
 * \f]
 * which are described completely by their waiting time distribution. This means that the time, the
 * random process \f$X_1\f$ needs to complete is gamma-distributed with shape \f$k=10\f$ and scale
 * \f$\theta=100\f$. Analogously, the processes \f$X_2\f$ and \f$X_3\f$ are defined by their
 * waiting time distribution.
 *
 * From these basic building blocks, a more complex random process can be assembled by
 * combining these blocks. Using the example above, one may define a new process that is a chain of
 * the processes \f$X_1\f$-\f$X_3\f$. Such a simple chain can be used to describe the successive
 * processing of information entering at time \f$t=0\f$ at the first stage described by \f$X_1\f$.
 * Once the first processing stage finished, its result gets forwarded to the second stage described
 * by \f$X_2\f$ and finally to the last stage described by \f$X_3\f$. The waiting time of the
 * complete chain is again a random variable that is the sum of all random variables,
 * i.e. mathematically
 * \f[
 *  Y = X_1 + X_2 + X_3\,.
 * \f]
 *
 * SochBB determines the probability density function (PDF) or cummulative probability function
 * (CDF) of the waiting-time distribution of the process \f$Y\f$ analytically (as far as
 * possible) or numerically (where an analytic method fails). More over it provides an efficient
 * and correct sampler for the random variables.
 *
 * As the summands commute and the distribution of the
 * sum of \f$X_1\f$ and \f$X_3\f$ can be determined analytically as
 * \f$X' = (X_1+X_3)\sim\Gamma(11,100)\f$, the random variable \f$Y\f$ can now be expressed as
 * \f$Y = X_2 + X'\f$ and only a single numeric convolution is necessary to obtain the PDF of the
 * random variable \f$Y\f$.
 *
 *
 * \section cli Command line interface
 * The StochBB installation comes with a command line tool. With this tool, a process definition
 * in @ref xml "XML" can be analyzed. The results are returned as CSV or they can be plotted
 * directly.
 *
 * \subsection clisyn Synopsis
 * \code
 *  stochbb [OPTIONS] INPUTFILE [OUTPUT OPTIONS]
 * \endcode
 *
 * \subsubsection cliops Options
 *  - `--help` Prints a short help string and exits.
 *  - `--version` Prints the version string and exits.
 *  - `--log-debug` Prints debug messages to stderr. By default, only warnings and error messages
 *    are printed to stderr.
 *
 * \subsubsection clioops Output options
 *  - `--pdf` Specifies to evaluate the PDF of the variables selected in the `INPUTFILE` (default).
 *  - `--cdf` Specifies to evaluate the CDF instead of the PDF of the variables selected in the
 *    `INPUTFILE`.
 *  - `--plot` Plots the PDF or CDF of the output variables specified in the file `INPUTFILE`.
 *  - `--csv=FILENAME` Writes the PDF or CDF of the output variables specified in the `INPUTFILE`
 *    to `FILENAME`.
 *
 * \subsection cliuse Usage
 *
 * \subsection cliex Examples
 *
 */

#ifndef __SBB_STOCHBB_HH__
#define __SBB_STOCHBB_HH__

#include "api.hh"

#endif // __SBB_STOCHBB_HH__

