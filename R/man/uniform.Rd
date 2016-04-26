\name{uniform}
\alias{uniform}
\title{
A uniformly distributed random variable
}
\description{
Constructs a uniformly distributed random variable.
}
\usage{
 uniform(a, b, name="")
}
\arguments{
 \item{a}{the lower bound.}
 \item{b}{the upper bound.}
 \item{name}{optional name of the random variable.}
}
\seealso{
 \code{\link{atomic}}
}
\examples{
\dontrun{

# Random variable X ~ U[0,1]
X <- uniform(0,1)
}
}
