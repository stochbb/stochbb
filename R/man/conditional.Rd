\name{conditional}
\alias{conditional}
\title{
Conditional random variables.
}
\description{
This function constructs a new random variable that is a mixture of two independent random variables
\eqn{Y_1, Y_2}{Y1,Y2}, where the relative weight of the random variables is determined by the
probability of \eqn{X_1 < X_2}{X1 < X2}. That is

\deqn{Z = Y_1\;if\,X_1 < X_2;\;Y_2\,else.}{Z = Y1 if X1<X2; Y2 else.}
}
\usage{
 conditional(X1, X2, Y1, Y2)
}
\arguments{
 \item{X1, X2}{the random variables to condition on.}
 \item{Y1, Y2}{the random variables to select from.}
}
\seealso{
 \code{\link{condchain}},
 \code{\link{independent}},
 \code{\link{derived}}.
}
