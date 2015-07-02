\name{logSum}
\alias{logSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate log(x + y) from log(x) and log(y)}
\description{
Given \eqn{\log(x)}{log(x)} and \eqn{\log(y)}{log(y)} this function calculates \eqn{\log(x + y)}{log(x + y)}
using the identity
\deqn{\log(x+y) = \log(x) + \log\left(1 + e^{\log(y) - \log(x)}\right)}{#
log(x + y) = log(x) + log(1 + exp[log(y) - log(x)])}
}
\usage{
logSum(x, y = NULL, base = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Numeric vector or matrix.}
  \item{y}{ Numeric vector or matrix of same dimensions as \code{x} or missing. }
  \item{base}{Base of the logarithm to be used. If \code{base} is 0 the natural logarithm is used.}
}
\details{
  If \code{y} is missing the function is applied recuresively to all elements of \code{x}, i.e., all
  elements of \code{x} are added up.  If both \code{x} and \code{y} are given they are added element wise.

}
\value{
  If only \code{x} is given a scalar is returned, representing the sum of all elements of \code{x}.
  Otherwise a vector or matrix of the same dimensions as \code{x} and \code{y}.
}
\author{Peter Humburg}
\note{This function is useful for cases where \code{x} and \code{y} cannot be represented accurately by machine numbers
but log(x) and log(y) can.}

\examples{
    x <- 1:4
    y <- 4:1
    ## calculate sum of x an y directly
    xy.s <- x + y
    ## and after log transformation
    xy.ls <- logSum(log(x),log(y))
    ## errors are small:
    err <- xy.ls - log(xy.s) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
