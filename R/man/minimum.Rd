\name{minimum}
\alias{minimum}
\title{
Represents the minimum of random independent variables.
}
\description{
This function constructs a new random variable that is the minimum of the given independent random
variables. That is \eqn{Z = \min(X_1,...,X_N)}{Z = min(X_1,...,X_N)}.
}
\usage{
 minimum(X1, X2, ...)
}
\arguments{
 \item{X1, X2, ...}{the random variables to get the minimum from.}
}
\seealso{
 \code{\link{maximum}}
 \code{\link{derived}}
}
\examples{
\dontrun{
# Z = min(X1, X2)
Z <- minimum(X1, X2)

}
}
