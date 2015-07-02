\name{viterbi}
\alias{viterbi}
\alias{viterbi,hmm-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Most Likely State Sequence Using the Viterbi Algorithm}
\description{
 The Viterbi algorithm computes the most likely sequence of states given
 an HMM and an observation sequence.
}
\usage{
\S4method{viterbi}{hmm}(hmm, obs, names=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm}{ Object of class \code{hmm}. }
  \item{obs}{ A vector containing the observation sequence. }
  \item{names}{ Logical indicating whether state names should be returned. 
  If \code{TRUE} (the default) the returned sequence consists of state names, 
  otherwise the state index is returned instead.}
}

\value{
  A list with components
  \item{stateSeq}{Most likely state sequence.}
  \item{logProb}{The probability of \code{stateSeq} given \code{hmm} and \code{obs}.}
  \item{matrix}{The dynamic programming matrix.}
}
\references{
Viterbi, A. J. 1967
Error bounds for convolutional codes and an assymptotically optimal decoding algorithm.
\emph{IEEE Transactions on Information Theory}, \bold{13}, 2600--269.
}
\author{Peter Humburg}

\seealso{\code{\link{viterbiTraining}}}
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
vit.res <- viterbi(model, obs$observation)

## how well did we do?
sum(vit.res$stateSeq == obs$states)/length(vit.res$stateSeq)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}

