\name{invgamma}
\alias{invgamma}
\title{
A inverse gamma distributed random variable
}
\description{
Constructs a inverse gamma distributed random variable. Either as an atomic random variable with fixed
parameters \code{alpha} & \code{beta}, or as a compound random variable if both or one of the parameters is a random
variable too.
}
\usage{
 invgamma(alpha, beta, name="")
}
\arguments{
 \item{alpha}{the shape parameter.}
 \item{beta}{the scale parameter.}
 \item{name}{optional name of the random variable.}
}
\seealso{
 \code{\link{atomic}}
}
\examples{
\dontrun{

# Random variable X ~ InvGama(1,1)
X <- invgamma(1,1)

# Compound random variable Y ~ invgamma(1,X)
Y <- invgamma(1,X)
}
}
