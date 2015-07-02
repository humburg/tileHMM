\name{forward}
\alias{forward}
\alias{backward}
\alias{forward,hmm-method}
\alias{backward,hmm-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of Forward and Backward Variables}
\description{
  These functions calculate the forward and backward variables for a given model
  and observation sequence.  All computations are carried out in log-space.
}
\usage{
\S4method{forward}{hmm}(hmm, obs)
\S4method{backward}{hmm}(hmm, obs)
}

\arguments{
  \item{hmm}{An object of class \code{hmm} or one of its subclasses representing the hidden Markov model.}
  \item{obs}{A vector containing the observation sequence.}
}

\value{
    \code{backward} returns the \eqn{N \times T} matrix of (log transformed) backward variables,
	where \eqn{N} is the number of states of \code{hmm} and \eqn{T} is the length of \code{obs}. 
    
    \code{forward} returns a list with components
  \item{logProb}{log[P(\code{obs}|\code{hmm})]}
  \item{alpha.scaled}{The matrix of log transformed forward variables. This has the same dimensions
  as the matrix returned by \code{backward}}
}
\references{
Rabiner, L. R. 1989
A tutorial on hidden Markov models and selected applications in speech recognition.
\emph{Proceedings of the IEEE}, \bold{77}(2), 257--286.
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
obs <- sampleSeq(model, 100)

## calculate the probability of the observation given the model
fwd <- forward(model, obs)
fwd$logProb

## compute posterior probabilities
bwd <- backward(model,obs)
post <- bwd + fwd$alpha.scaled
post <- t(t(post) - apply(post,2,logSum))

## get sequence of most likely states
state.seq <- state.names[apply(post,2,which.max)]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}