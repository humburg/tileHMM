\name{viterbiTraining}
\alias{viterbiTraining}
\alias{viterbiTraining,hmm,list-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate HMM Parameters Using Viterbi Training}
\description{
  Viterbi training is a faster but less reliable alternative to Baum-Welch for parameter estimation.
}
\usage{
\S4method{viterbiTraining}{hmm,list}(hmm, obs, max.iter=10, eps=0.01, 
        df=NULL, trans.prior=NULL, init.prior=NULL, keep.models=NULL, verbose=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm}{ Object of class \code{hmm}. }
  \item{obs}{ List of observation sequences. }
  \item{max.iter}{ Maximum number of iterations. }
  \item{eps}{ Minimum change in log likelihood between successive iterations.}
  \item{df}{If this is \code{NULL} the degrees of freedom for the t distributions are e
  stimated from the data.  Otherwise they are set to \code{df}.}
  \item{trans.prior}{Prior distribution of transition probabilities. A prior can be specified
  either by providing a matrix with transition probabilities or by setting \code{trans.prior=TRUE}. 
  In the latter case the initial parameter estimates are used as prior. If \code{trans.prior} is \code{NULL}
  (the default) no prior is used.}
  \item{init.prior}{Prior distribution of initial state probabilities. A prior can be specified
  either by providing a vector with initial state probabilities or by setting \code{init.prior=TRUE}. 
  In the latter case the initial parameter estimates are used as prior. If \code{init.prior} is \code{NULL}
  (the default) no prior is used.}
  \item{keep.models}{A character string interpreted as a file name. If \code{keep.models} is not \code{NULL}
  the models produced during the parameter estimation procedure are saved to a file.}
  \item{verbose}{Level of verbosity. Allows some control over the amount of output 
  printed to the console.}
}
\value{
  Object of class \code{hmm} with the pest parameter estimates (in terms of likelihood) found
  during the fitting procedure.
}
\references{
Juang, B.-H. and Rabiner, L. R. 1990
A segmental k-means algorithm for estimating parameters of hidden Markov models. 
\emph{IEEE Transactions on Acoustics, Speech, and Signal Processing}, \bold{38}(9), 1639--1641.
}
\author{Peter Humburg}

\seealso{\code{\link{viterbi}}, \code{\link{baumWelch}}, \code{\link{viterbiEM}}, \code{\link{hmm.setup}}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.035, 0.01)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
hmm1 <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df), 
    state.names)

## generate observation sequences from model
obs.lst <- list()
for(i in 1:50) obs.lst[[i]] <- sampleSeq(hmm1, 100)

## fit an HMM to the data (with fixed degrees of freedom)
hmm2 <- hmm.setup(obs.lst, state=c("one","two"), df=5)
hmm2.fit <- viterbiTraining(hmm2, obs.lst, max.iter=20, df=5, verbose=1)

## fit an HMM to the data, this time estimating the degrees of freedom
hmm3.fit <- viterbiTraining(hmm2, obs.lst, max.iter=20, verbose=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
\keyword{models}
