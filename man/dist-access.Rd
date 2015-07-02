\name{dist-access}
\alias{as.data.frame.discDist}
\alias{as.vector.discDist}
\alias{as.matrix.discDist}
\alias{length,discDist-method}
\alias{[.discDist}
\alias{[[.discDist}
\alias{[.contDist}
\alias{[<-.contDist}
\alias{[.tDist}
\alias{[<-.tDist}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessing and Converting Objects of Class "dist"}
\description{
These methods provide convenient access to objects of class \code{dist} as well as 
conversion to other data structures.
}
\usage{
\method{as.data.frame}{discDist}(x, row.names = NULL, optional = FALSE, ...)
\method{as.vector}{discDist}(x, mode = "any")
\method{as.matrix}{discDist}(x, ...)

\S4method{length}{discDist}(x)

\method{[}{discDist}(x, i, ...)
\method{[[}{discDist}(x, i, ...)
\method{[}{contDist}(x, i, j, ...)
\method{[}{tDist}(x, i, j, log=FALSE, ...)

\method{[}{contDist}(x, i, j) <- value
\method{[}{tDist}(x, i, j) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{discDist}, \code{contDist} or \code{tDist}.}
  \item{i}{If \code{x} is an object of class \code{discDist} this is expected to be
      a symbol from the alphabet of \code{x}. For all other classes \code{i} is interpreted
      as the index of the mixture component that should be accessed. \code{i} may be missing
      in which case values for all mixture components are returned or replaced.}
  \item{j}{Either a character string identifying one of the parameters of \code{x} or a numeric value.
      In the later case the density at point \code{j} is returned.}
  \item{log}{Logical indicating whether the density value should be log transformed.} 
  \item{row.names}{NULL or a character vector giving the row names for the data frame. Missing values are not allowed.}
  \item{optional}{logical. If TRUE, setting row names and converting column names (to syntactic names) is optional.}
  \item{mode}{Mode of \code{vector}.}
  \item{value}{New parameter value(s).}
  \item{\dots}{Additional arguments.}
}

\value{
  The conversion functions return objects of the respective class.
  
  \code{length} returns the number of symbols in the alphabet of \code{x}.
  
  The access functions return either the requested parameter value or the value of the dansity
  function at point \code{j}.
}
\author{Peter Humburg}

\seealso{\code{\linkS4class{dist}}}
\examples{
## converting "discDist" objects
nucleo <- new("discDist",alpha=c("A","C","G","T"),prob=rep(0.25,times=4))
nucleo.vec <- as.vector(nucleo)
nucleo.mat <- as.matrix(nucleo)
nucleo.df <- as.data.frame(nucleo)

## get number of symbols in alphabet
len <- length(nucleo)   # = 4

## get probability for symbol 'A'
nucleo['A']  # = 0.25

## accessing tDist objects
mydt <- new("tDist", mean=0, var=1, df=3)

## evaluate density function
mydt[, 2]  # = 0.06750966
mydt[, 1000] # = 3.307954e-12
mydt[, 1000, log=TRUE]  # = -26.43469

## access parameter values
mydt[, "mean"]   # = 0

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

