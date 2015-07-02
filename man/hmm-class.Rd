\name{hmm-class}
\docType{class}
\alias{hmm-class}
\alias{show,hmm-method}

\title{Class "hmm"}
\description{Virtual base class for HMMs.}
\section{Objects from the Class}{
Do not create objects of this class directly.  Instead use objects of derived classes like
\code{"contHMM"}.
}
\section{Slots}{
	 \describe{
    \item{\code{transition.matrix}:}{Object of class \code{"matrix"}, storing the transition probabilities of the Markov chain.}
    \item{\code{emission}:}{Object of class \code{"list"} containing objects of class \code{"dist"} to represent emission distributions for each state.}
    \item{\code{init}:}{Object of class \code{"numeric"}. The initial state distribution of the Markov chain.}
  }
}
\section{Methods}{
  \describe{
    \item{baumWelch}{\code{signature(hmm = "hmm", obs = "list")}: Baum-Welch algorithm to estimate model parameters.}
    \item{viterbiTraining}{\code{signature(hmm = "hmm", obs = "list")}: Viterbi training to estimate model parameters.}
    \item{viterbi}{\code{signature(hmm = "hmm")}: Viterbi algorithm to calculate most likely state sequence.}
    \item{backward}{\code{signature(hmm = "hmm")}: Computing of backward variables.}
    \item{forward}{\code{signature(hmm = "hmm")}: Computing of forward variables.}
    \item{states}{\code{signature(hmm = "hmm")}: Returns state names.}
    \item{show}{\code{signature(object = "hmm")}: Print summary of hmm object.}
	 }
}
\author{ Peter Humburg }

\seealso{
	\code{\linkS4class{contHMM}}, \code{\link{baumWelch}}, \code{\link{viterbiTraining}},
	\code{\link{viterbi}}, \code{\link{forward}}, \code{\link{backward}}, \code{\link{states}}
}
\examples{
showClass("hmm")
}
\keyword{classes}
\keyword{models}
