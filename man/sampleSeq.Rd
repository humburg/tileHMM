\name{sampleSeq}
\alias{sampleSeq}
\alias{sampleSeq,contHMM,numeric-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate Observation Sequence from HMM}
\description{
 Generates an observation sequence according to a model given as object
 of class \code{hmm}. Optionally the undelying state sequence is returned
 together with the observations.
}
\usage{
\S4method{sampleSeq}{contHMM,numeric}(hmm, size, return.states=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm}{Object of class \code{contHMM}.}
  \item{size}{Numeric value indicating the desired length of the observation sequence.}
  \item{return.states}{Logical indicating whether the underlying state sequence should be returned
  together with the observation sequence.}
}

\value{
  If \code{return.states} is \code{FALSE} (the default) a vector of length \code{size} with observations
  sampled from \code{hmm}. If \code{return.states} is \code{TRUE} a list with components
  \item{states}{Character vector with state sequence that was used to generate the observations.}
  \item{observation}{Vector with observations sampled from \code{hmm}.}
}
\author{Peter Humburg}
\seealso{\code{\linkS4class{hmm}}, \code{\link{states}}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.035, 0.01)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
hmm1 <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df),
    state.names)

## generate observation sequence from model
obs <- sampleSeq(hmm1, 100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}