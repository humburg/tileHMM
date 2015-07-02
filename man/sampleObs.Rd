\name{sampleObs}
\alias{sampleObs}
\alias{sampleObs,discDist,numeric-method}
\alias{sampleObs,tDist,numeric-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample Observations from Probability Distribution}
\description{
  Draws a sample of given size from a discrete or continuous probability distribution.
}
\usage{
\S4method{sampleObs}{discDist,numeric}(dist, size)
\S4method{sampleObs}{tDist,numeric}(dist, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{Object of class \code{discDist} or \code{tDist}. The probability distribution to use.}
  \item{size}{Sample size.}
  \item{\dots}{Additional arguments.}
}

\value{
  Vector of sampled values.
}
\author{Peter Humburg}

\seealso{\code{\linkS4class{dist}}}
\examples{
## sampling from a t distribution
tdist <- new("tDist",mean=0,var=1,df=3)
obs <- sampleObs(tdist,100)

## sampling from a discrete distribution
nucleo <- new("discDist",alpha=c("A","C","G","T"),prob=rep(0.25,times=4))
dna <- sampleObs(nucleo,100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
