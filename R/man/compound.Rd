\name{compund}
\alias{compound}
\title{
Compound random variables
}
\description{
Compound random variables are \link{derived} random variables are distributed according to some
parametric density function (e.g., \code{\link{normal}}, \code{\link{gamma}}, etc.) where at least 
one of the distribution parameters is a random variable too. In \pkg{stochbb} compound random 
variables are constructed using the same factory functions as the \link{atomic} random variables 
(see below) but passing a random variable for at least one of the distribution parameters. 

}
\seealso{
 \code{\link{delta}},
 \code{\link{uniform}},
 \code{\link{normal}},
 \code{\link{gamma}},
 \code{\link{invgamma}},
 \code{\link{weibull}},
 \code{\link{derived}}.
}
