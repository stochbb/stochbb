\name{studt}
\alias{studt}
\title{
A Student't distributed random variable
}
\description{
Constructs a Student's t distributed random variable. Either as an atomic random variable with fixed
parameters nu, or as a compound random variable if the the parameter is a random
variable too.
}
\usage{
 studt(nu, name="")
}
\arguments{
 \item{nu}{the degrees of freedom.}
 \item{name}{optional name of the random variable.}
}
\seealso{
 \code{\link{atomic}}
}
\examples{
\dontrun{

# Random variable X ~ t(1) (a.k.a cauchy)
X <- studt(1)
}
}
