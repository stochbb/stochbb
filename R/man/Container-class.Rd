\name{Container-class}
\Rdversion{1.1}
\docType{class}
\alias{Container-class}
\alias{Container}
\alias{Rcpp_Container}
\title{Base for all \pkg{stochbb} classes}
\description{
This is simply the base class for all classes in \pkg{stochbb}.
}
\section{Extends}{
Class \code{\linkS4class{C++Object}}, directly.

}
\seealso{
 \code{\linkS4class{Var}},
 \code{\linkS4class{Density}},
 \code{\linkS4class{ExactSampler}}.
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{isNull}:}{Returns TRUE if the container is empty.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{repr()}:}{ Returns a string representation of the held object.}
  }
}
