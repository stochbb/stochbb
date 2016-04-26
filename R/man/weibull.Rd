\name{weibull}
\alias{weibull}
\title{
A Weibull distributed random variable
}
\description{
Constructs a Weibull distributed random variable. Either as an atomic random variable with fixed
parameters k & lambda, or as a compound random variable if both or one of the parameters is a random
variable too.
}
\usage{
 weibull(k, lambda, name="")
}
\arguments{
 \item{k}{the shape parameter.}
 \item{lambda}{the scale parameter.}
 \item{name}{optional name of the random variable.}
}
\seealso{
 \code{\link{atomic}}
}
\examples{
\dontrun{

# Random variable X ~ Weibull(1,1)
X <- weibull(1,1)

# Compound random variable Y ~ Webull(X,1)
Y <- weibull(1,X)
}
}
