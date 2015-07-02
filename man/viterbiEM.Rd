\name{viterbiEM}
\alias{viterbiEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Efficient Estimation of HMM Parameters}
\description{
  Uses a combination of Viterbi training and Baum-Welch algorithm to estimate parameters for a hidden Markov model.
}
\usage{
viterbiEM(hmm, data, max.iter = c(5, 15), eps = 0.01, 
    verbose = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm}{Object of class \code{hmm}. This is used as starting point for the optimisation procedure.}
  \item{data}{A \code{list} of observation sequences.}
  \item{max.iter}{Maximum number of iterations (see Details).}
  \item{eps}{Minimum change in log-likelihood between iterations (see Details). }
  \item{...}{Additional arguments to be passed to \code{viterbiTraining} and \code{baumWelch} (see Details).}
  \item{verbose}{Level of verbosity. Higher numbers produce more status messages.}
}
\details{
  The values of arguments \code{max.iter} and \code{eps} can have either one or two elements.  In the latter case
  the first element is used for \code{viterbiTraining} and the second one for \code{baumWelch}.
  
  Additional arguments can be passed to \code{viterbiTraining} and \code{baumWelch} by using arguments of the form
  \code{viterbi = list(a = a.value)} and \code{baumWelch = list(b = b.value)} respectively. All other arguments are 
  passed on to both functions. 
}
\value{
 An object of class \code{hmm} with optimised parameter estimates.
}
\references{ 
Humburg, P. and Bulger, D. and Stone, G.
Parameter estimation for robust HMM analysis of ChIP-chip data. BMC Bioinformatics 2008, 9:343}
\author{Peter Humburg}

\seealso{\code{\link{baumWelch}}, \code{\link{viterbiTraining}}, \code{\link{hmm.setup}}}
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
hmm2.fit <- viterbiEM(hmm2, obs.lst, max.iter=c(5,15), verbose=2, df=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
