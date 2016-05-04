\name{Density-class}
\Rdversion{1.1}
\docType{class}
\alias{Density-class}
\alias{Density}
\alias{Rcpp_Density}
\title{Represents a density of a random variable.}
\description{
Instances of this class are returned by the \code{density} field of \code{\linkS4class{Var}}
instances and represents the density of the random variable.

Please note that there is no direct constructor for densities. Instances of this class are always
associated with an instance of a random variable.
}
\section{Extends}{
Class \code{\linkS4class{Container}}, directly.

}
\seealso{
\code{\linkS4class{Var}}
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{isNull}:}{Is TRUE if the density is invalid.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{repr()}:}{Returns a string representation of the object.}
    \item{\code{rangeEst(alpha)}:}{Returns a rough estimate of the \code{alpha}-quantiles of the
    density. This can be used to determine the interval on which the density should be evaluated.}
    \item{\code{evalCDF(Tmin, Tmax, out)}:}{Evaluates the cumulative distribution function
    (probability function) of the density on the interval \eqn{[T_{min},T_{max}]}{[Tmin,Tmax]} and
    stores the result in the vector \code{out}. The dimension of the vector specifies the number
    of bins for the evaluation of the CDF.}
    \item{\code{eval(Tmin, Tmax, out)}:}{Evaluates the probability density function
    (PDF) on the interval \eqn{[T_{min},T_{max}]}{[Tmin,Tmax]} and stores the result in the vector
    \code{out}. The dimension of the vector specifies the number of bins for the evaluation of the
    PDF.}
  }
}
