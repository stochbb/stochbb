\name{affine}
\alias{affine}
\title{
Affine transformation of a random variable.
}
\description{
This function constructs a new random variable that is the affine transformed of the given random
variable. That is \eqn{Y = a\,X+b}{Y = a*X + b}.

Please note that the \code{\%+\%} operator can also be used to shift a random variable. That is,
\code{Z <- X \%+\% 100} is equivalent to \code{Z <- affine(X, 1, 100)}.
}
\usage{
 affine(X, a, b)
}
\arguments{
 \item{X}{the random variable.}
 \item{a}{the scaleing.}
 \item{b}{the shift.}
}
\seealso{
 \code{\link{derived}}
}
\examples{
\dontrun{

# Z = 3*X + 100
Z <- affine(X, 3, 100)

# Z = X + 100
Z <- affine(X, 1, 100)
# or equivalently
Z <- X \%+\% 100
}
}
