\name{Var-class}
\Rdversion{1.1}
\docType{class}
\alias{Var-class}
\alias{Var}
\alias{Rcpp_Var}
\alias{is.Var}
\title{Class \code{Var} for all random variables}
\description{
Represents a random variable within a system of dependent random variables. Please note that there
is no direct constructor for the \code{Var} class. Instances of this class are created by
construing an \link{atomic} or \link{derived} random variable using factory functions.
}
\section{Extends}{
Class \code{\linkS4class{Container}}, directly.
}

\seealso{
 \code{\linkS4class{Density}}
 \code{\linkS4class{Container}}
 \code{\link{atomic}}
 \code{\link{derived}}
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{density}:}{Holds the density of the random variable.}
    \item{\code{isNull}:}{If \code{TRUE}, the container is empty.}
    \item{\code{name}:}{Holds the name of the random variable.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{repr()}:}{ Returns a string representation of the random variable. }
    \item{\code{mutuallyIndep(x)}:}{ Returns true if the given random variable is independent from
    the this one.}
    \item{\code{dependsOn(x)}:}{ Returns true if this random variable depends on the given one.}
  }
}


