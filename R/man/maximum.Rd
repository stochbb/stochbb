\name{maximum}
\alias{maximum}
\title{
Represents the maximum of random independent variables.
}
\description{
This function constructs a new random variable that is the maximum of the given independent random
variables. That is \eqn{Z = \max(X_1,...,X_N)}{Z = max(X_1,...,X_N)}.
}
\usage{
 maximum(X1, X2, ...)
}
\arguments{
 \item{X1, X2, ...}{the random variables to get the maximum from.}
}
\seealso{
 \code{\link{minimum}},
 \code{\link{derived}}.
}
\examples{
\dontrun{

# Z = max(X1, X2)
Z <- maximum(X1, X2)
}
}
