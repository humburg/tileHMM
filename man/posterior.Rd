\name{posterior}
\alias{posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Posterior Probability for States of HMM}
\description{
  For each state of an HMM the posterior probability that this state produced a given observation is calculated.
}
\usage{
posterior(data, hmm, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector with observation sequence.}
  \item{hmm}{Object of class \code{hmm}.}
  \item{log}{Logical indicating whether the logarithm of the posterior probability should be returned.}
}
\details{
  Regardless of the value of \code{log} the computation is carried out in log space. If \code{log = FALSE}
  the result is transformed back to linear space before it is returned.  
}
\value{
  A matrix with as many rows as \code{hmm} has states and one column for each entry in data.
}
\author{Peter Humburg}

\seealso{\code{\link{forward}}, \code{\link{backward}}}
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

## calculate posterior probability for state "one"
post <- posterior(obs, model, log=FALSE)

## get sequence of individually most likely states
state.seq <- apply(post,2,max)
state.seq <- states(model)[state.seq]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
