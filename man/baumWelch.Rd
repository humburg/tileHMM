\name{baumWelch}
\alias{baumWelch}
\alias{baumWelch,hmm,list-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Baum-Welch Algorithm}
\description{
The Baum-Welch algorithm [Baum et al., 1970] is a well established method for 
estimating parameters of HMMs.  It represents the EM algorithm [Dempster et al., 1977] 
for the specific case of HMMs.  The formulation of the Baum-Welch algorithm 
used in this implementation is based on the description given by Rabiner [1989].
}
\usage{
\S4method{baumWelch}{hmm,list}(hmm, obs, max.iter=FALSE, eps=0.01, df=NULL, 
	trans.prior=NULL, init.prior=NULL, verbose=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm}{An object of class \code{hmm} or one of its subclasses representing the hidden Markov model.}
  \item{obs}{A list of observation sequences.}
  \item{max.iter}{Maximum number of iterations. (optional)}
  \item{eps}{Minimum difference in log likelihood between iterations. Default: 0.01}
  \item{df}{If this is \code{NULL} the degrees of freedom for the t distributions are estimated 
  from the data.  Otherwise they are set to \code{df}.}
  \item{trans.prior}{Prior distribution of transition probabilities. A prior can be specified
  either by providing a matrix with transition probabilities or by setting \code{trans.prior=TRUE}. 
  In the latter case the initial parameter estimates are used as prior. If \code{trans.prior} is \code{NULL}
  (the default) no prior is used.}
  \item{init.prior}{Prior distribution of initial state probabilities. A prior can be specified
  either by providing a vector with initial state probabilities or by setting \code{init.prior=TRUE}. 
  In the latter case the initial parameter estimates are used as prior. If \code{init.prior} is \code{NULL}
  (the default) no prior is used.}
  \item{verbose}{Level of verbosity. Allows some control over the amount of output 
  printed to the console.}
}

\value{
Returns the HMM with optimised parameters.
}
\references{
Baum, L. E. and Petrie, T. and Soules, G. and Weiss, N. 1970 
A maximization technique occuring in the statistical analysis of 
probabilistic functions of markov chains. 
\emph{The Annals of Mathematical Statistics}, \bold{41}(1), 164--171.

Dempster, A. P. and Laird, N. M. and Rubin, D. B. 1977 
Maximum likelihood for incomplete data via the EM algorithm. 
\emph{Journal of the Royal Statistical Society}, Series B, \bold{39}(1).

Rabiner, L. R. 1989
A tutorial on hidden Markov models and selected applications in speech recognition.
\emph{Proceedings of the IEEE}, \bold{77}(2), 257--286.
}
\author{Peter Humburg}

\seealso{\code{\link{viterbiTraining}}, \code{\link{viterbiEM}}, \code{\link{getHMM}}, \code{\link{hmm.setup}} }
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
hmm2.fit <- baumWelch(hmm2, obs.lst, max.iter=20, df=5, verbose=1)

## fit an HMM to the data, this time estimating the degrees of freedom
hmm3.fit <- baumWelch(hmm2, obs.lst, max.iter=20, verbose=1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
\keyword{models}
