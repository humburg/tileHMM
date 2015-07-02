\name{states}
\alias{states}
\alias{states,hmm-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{State Names of Hidden Markov Model}
\description{
Returns a vector of state names for given HMM.
}
\usage{
\S4method{states}{hmm}(hmm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm}{Object of class \code{hmm}}
  \item{\dots}{Additional arguments.}
}

\value{
A character vector with as many entries as \code{hmm} has states. Each entry represents the name of the
corresponding state.
}
\author{Peter Humburg}

\seealso{\code{\linkS4class{hmm}}}
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
obs <- sampleSeq(model, 100, return.states=TRUE)

## compute most likely state sequence for obs
## return sequence of state indices instead of names
vit.res <- viterbi(model, obs$observation, names=FALSE)

## get sequence of state names via call to 'states'
state.seq <- states(model)[vit.res$stateSeq]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
