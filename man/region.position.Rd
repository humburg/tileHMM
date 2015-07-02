\name{region.position}
\alias{region.position}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify Enriched Regions}
\description{
After calling individual probes enriched or non-enriched this function can be used
to combine probes into enriched regions. 
}
\usage{
region.position(probe.calls, region = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probe.calls}{A vector of probe calls.}
  \item{region}{Probe call that should be combined into regions.}
}

\value{
A matrix with two rows and as many columns as there are enriched regions.
The first row gives the index of the first probe of each region the second 
row the index of the last probe inside each region.
}
\author{Peter Humburg}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.1, 0.02)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
model <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df), 
    state.names)

## obtain observation sequence from model
obs <- sampleSeq(model, 100)

## compute most likely state sequence for obs
vit.res <- viterbi(model, obs)

## find regions attributed to state "one"
reg.pos <- region.position(vit.res$stateSeq, region="one")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}