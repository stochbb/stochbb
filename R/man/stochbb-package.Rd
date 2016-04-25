\name{stochbb-package}
\alias{stochbb-package}
\alias{stochbb}
\docType{package}
\title{
\packageTitle{stochbb}
}
\description{
\packageDescription{stochbb}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{stochbb}
\packageIndices{stochbb}
This package provides a collection of classes and functions to assemble and analyze complex systems
of dependent random variables. To this end, it allows to analyze simple stochastic "processes" by
means of numerical approximation rather than stochastic simulation.

Many computational models of cognitive processes are described as "stochastic processes", although
they are formally not stochastic processes (neither in time nor space). Mathematically, they are
"just" a system of dependent random variables, whose distributions do not depend on time or space.

Many of these models are described as a system of "stages", each having a specified waiting time
distribution which might be modulated by the presented simulus or the processing time needed for a
previous stimulus. These stages are then chained, that is, as one stage completes it may trigger
some other stage. Thus forming a simple linear processing \code{\link{chain}} of some initial
trigger (stimulus) upto the final response stage.

Processing may also be performed in parallel, that is one stage triggers several other stages
which must then be combined to a single response: for example, waiting until all parallel stages are
finished (\code{\link{maximum}}) or until the first stage finishes (\code{\link{minimum}}). Parallel
stages may also inhibit each other: The completion of one stage may cancel the ongoing processing
along all other parallel paths (\code{\link{condchain}}).

The "gold standard" approach to analyze such systems of dependent random variables is stochastic
simulation. Although this is a simple and exact approach to sample from these random variables, it
is also slow. A large number of samples must be drawn from all involved random variables to obtain
a reliable statistic about the final response variable.

Many of the afroementioned operations on random variables, however, can be described analytically.
For example, the probability function (CDF) of a random variable \eqn{Z}{Z} being the maximum of
two independent random variables \eqn{X}{X} and \eqn{Y}{X}, \eqn{Z=\max(X,Y)}{Z = max(X, Y)}, is
\eqn{F_Z(z)=F_X(x)\,F_Y(y)}{Fz(z) = Fx(x) * Fy(y)}, where \eqn{F_X}{Fx} and \eqn{F_Y}{Fy} are the
probability functions of the random variables \eqn{X}{X} and \eqn{Y}{Y}, respectively. Consequently,
the probability density function \eqn{f_Z}{fz} of \eqn{Z}{Z} is then
\eqn{f_Z(z)=f_X(x)\,F_Y(y)+F_X(x)\,f_y(y)}{Fz(z) = fx(x) * Fy(y) + Fx(x) * fy(y)}. Similar relations
can be obtained for the minimum of random variables, mixtures of random variables and even for some
special cases of chained stages (sums of independent random variables). For the latter case, at
least a numerical approximation of the PDF by means of FFT convolution can be obtained efficiently.

The \pkg{stochbb} package allows to infer the marginal distributions of a system of dependent random
variables analytically and resorts to efficient numerical approximations when necessary. To this
end, it allows to analyze such models with high precision without resorting to expensive stochastic
simulation procedures.

Of cause it also provides an exact sampler for verification, if desired (see
\code{\linkS4class{ExactSampler}}).
}
\author{
\packageAuthor{stochbb}

Maintainer: \packageMaintainer{stochbb}
}
\keyword{ package }
