\name{gamma}
\alias{gamma}
\title{
A gamma distributed random variable
}
\description{
Constructs a gamma distributed random variable. Either as an atomic random variable with fixed
parameters \code{k} & \code{theta}, or as a compound random variable if both or one of the 
parameters is a random variable too.
}
\usage{
 gamma(k, theta, name="")
}
\arguments{
 \item{k}{the shape parameter.}
 \item{theta}{the scale parameter.}
 \item{name}{optional name of the random variable.}
}
\seealso{
 \code{\link{atomic}}
}
\examples{
\dontrun{

# Random variable X ~ Gamma(1,1)
X <- gamma(1,1)

# Compound random variable Y ~ gamma(1,X)
Y <- gamma(1,X)
}
}
