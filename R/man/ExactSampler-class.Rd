\name{ExactSampler-class}
\Rdversion{1.1}
\docType{class}
\alias{ExactSampler-class}
\alias{ExactSampler}
\alias{Rcpp_ExactSampler}
\title{Implements an exact sampler.}
\description{
This class implements a reasonable fast exact sampler for one or more possibly dependent random
variables.

}
\section{Extends}{
Class \code{\linkS4class{Container}}, directly.

}

\seealso{
  \code{\linkS4class{Var}}
  \code{\link{independent}}
}
\examples{
\dontrun{

 X <- normal(0,1)
 Y <- X \%+\% normal(0,1)
 sam <- new(ExactSampler, c(X,Y))
 res <- array(0, c(1000, 2))
 sam$sample(res)
}
}

\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{isNull}:}{Is \code{TRUE} if the container is invalid.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{initialize(vars)}:}{Constructs a new sampler. \code{vars} is the list of
    random variables to sample from.}
    \item{\code{sample(out)}:}{Samples from the random variable passed to the constructor. The
    number of columns of \code{out} must match the number of variables passed to the constructor,
    while the number of rows of \code{code} specifies the number of samples to draw.}
    \item{\code{repr()}:}{Returns a string representation of the object.}
  }
}
