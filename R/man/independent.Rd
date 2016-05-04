\name{independent}
\alias{independent}
\title{
Test for independence.
}
\description{
This function returns \code{TRUE} if all given random variables are mutually independent. Please
note that all \link{atomic} random variables are always independent from each other.

}
\usage{
 independent(X1, X2, ...)
}
\arguments{
 \item{X1, X2, ...}{the random variables to test.}
}
\examples{
\dontrun{

# Z = X1 + X2
X <- normal(0,1)
Y <- normal(0,1)
Z <- X \%+\% Y
# TRUE
independent(X,Y)

# FALSE
independent(X,Z)
}
}
