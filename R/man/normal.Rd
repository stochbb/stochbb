\name{normal}
\alias{normal}
\title{
A normal distributed random variable
}
\description{
Constructs a normal distributed random variable. Either as an atomic random variable with fixed
parameters \code{mu} & \code{sigma}, or as a compound random variable if both or one of the parameters is a random
variable too.
}
\usage{
 normal(mu, sigma, name="")
}
\arguments{
 \item{mu}{the mean parameter.}
 \item{sigma}{the standard deviation parameter.}
 \item{name}{optional name of the random variable.}
}
\seealso{
 \code{\link{atomic}}
}
\examples{
\dontrun{

# Random variable X ~ Norm(0,1)
X <- normal(0,1)

# Compound random variable Y ~ Norm(X,1)
Y <- normal(X,1)
}
}
