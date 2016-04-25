\name{Var-class}
\Rdversion{1.1}
\docType{class}
\alias{Var-class}

\title{Class \code{Var} for all random variables}
\description{
Represents a random variable within a system of dependent random variables. Please note that there
is no direct constructor for the \code{Var} class. Instances of this class are created by
construcing an \link{atomic} or \link{derived} random variable using functions.
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
    \item{\code{isNull}:}{If true, the container is empty.}
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


