\name{mixture}
\alias{mixture}
\title{
Represents a mixture of two independent random variables.
}
\description{
This function constructs a new random variable that is the mixture of the given independent random
variables.
}
\usage{
 mixture(w1, X1, w2, X2)
}
\arguments{
 \item{X1, X2}{the random variables to mix.}
 \item{w1, w2}{the relative weights of the random variables. They do not need to sum up to one.}
}
\seealso{
 \code{\link{derived}}
}
