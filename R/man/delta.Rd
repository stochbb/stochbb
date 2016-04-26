\name{delta}
\alias{delta}
\title{
A \emph{delta} distributed random variable
}
\description{
Constructs a \emph{delta} distributed random variable. That is a random variable that has a fixed value.
Such random variables are only used if a constant value is needed even if a random variable is
expected by some function or object. Such random variables are usually not used directly but are
needed internally, for example for compound distributed random variables where some parameters are
fixed.
}
\usage{
 delta(a)
}
\arguments{
 \item{a}{the fixed value of the \emph{random variable}.}
}
\seealso{
 \code{\link{atomic}}
}
