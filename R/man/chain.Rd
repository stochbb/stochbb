\name{chain}
\alias{chain}
\alias{\%+\%}
\title{
Represents a sum of independent random variables.
}
\description{
This function constructs a new random variable that is the sum of the given independent random
variables.

Please note that the \code{\%+\%} operator can also be used to shift a random variable. That is,
\code{Z <- X \%+\% 100} is equivalent to \code{Z <- affine(X, 1, 100)}.
}
\usage{
 chain(X1, X2, ...)
}
\arguments{
 \item{X1, X2, ...}{the random variables to sum.}
}
\seealso{
 \code{\link{affine}},
 \code{\link{independent}},
 \code{\link{derived}}.
}
\examples{
\dontrun{

# Z = X1 + X2
Z <- chain(X1, X2)

# or equivalently
Z <- X1 \%+\% X2
}
}
