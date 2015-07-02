\name{region.length}
\alias{region.length}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine Length of Positive and Negative Regions}
\description{
Given a logical vector indicating positive and negative probes
this function returns a list with components 'positive' and 'negative' providing length
information for positive and negative regions.
}
\usage{
region.length(probes, min.len=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{A logical vector indicating the position of enriched and non-enriched probes with 
  \code{TRUE} and \code{FALSE} respectively.}
  \item{min.len}{The minimal number of consecutive probes required to form a region.}
}

\value{
  A list with components \code{positive} and \code{negative}, each containing a numeric vector
  with the length of identified regions. 
}
\author{Peter Humburg}

\seealso{\code{\link{region.position}}}
\examples{
## create random probe calls
probes <- sample(c(TRUE,FALSE), 200, replace=TRUE)

## find length of all regions that contain at least two probes
reg.len <- region.length(probes, min.len=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
