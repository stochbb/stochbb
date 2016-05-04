\name{condchain}
\alias{condchain}
\title{
Chained conditional random variables.
}
\description{
This function constructs a new random variable similar to \code{\link{conditional}} but also adds
the minimum of the conditional random variables. That is

\deqn{Z = X_1+Y_1\;if\,X_1 < X_2;\;X_2+Y_2\,else.}{Z = X1+Y1 if X1<X2; X2+Y2 else.}
}
\usage{
 condchain(X1, X2, Y1, Y2)
}
\arguments{
 \item{X1, X2}{the random variables to condition on.}
 \item{Y1, Y2}{the random variables to select from.}
}
\seealso{
 \code{\link{conditional}}, 
 \code{\link{independent}},
 \code{\link{derived}}.
}
