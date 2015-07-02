\name{contHMM-class}
\docType{class}
\alias{contHMM-class}
\alias{initialize,contHMM}
%\alias{.baumWelchEmission,contHMM,contDist,list-method}
%\alias{.baumWelchEmission,contHMM,tDist,list-method}
%\alias{.baumWelchStep,contHMM,list-method}
%\alias{sampleSeq,contHMM,numeric-method}
%\alias{.viterbiTrainingEmission,contHMM,list,list-method}

\title{Class "contHMM"}
\description{Class for HMMs with continuous emission distributions.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("contHMM", transition, emission, init)}.
For the special case of \emph{t} distributions objects of this class can be created more
conveniently by a call to \code{getHMM}.  The function \code{hmm.setup} provides facilities to
create \code{"contHMM"} objects with initial parameter estimates obtained from data. 
}
\section{Slots}{
	 \describe{
    \item{\code{transition.matrix}:}{Object of class \code{"matrix"}, storing the transition probabilities of the Markov chain.}
    \item{\code{emission}:}{Object of class \code{"list"} containing objects of class \code{"contDist"} to represent emission distributions for each state.}
    \item{\code{init}:}{Object of class \code{"numeric"}. The initial state distribution of the Markov chain.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{hmm}"}, directly.
}
\section{Methods}{
  \describe{
   \item{sampleSeq}{\code{signature(hmm = "contHMM", size = "numeric")}: Sample observation sequence of length \code{size} from model }
	}
}

\author{Peter Humburg}

\seealso{
	\code{\linkS4class{contDist}}, \code{\linkS4class{tDist}},\code{\linkS4class{hmm}}, 
	\code{\link[=initialize,contHMM-method]{initialize}}, \code{\link{sampleSeq}}, 
	\code{\link{baumWelch}}, \code{\link{viterbiTraining}},	\code{\link{viterbi}}, 
	\code{\link{forward}}, \code{\link{backward}}, \code{\link{states}}
}
\examples{
showClass("contHMM")
}
\keyword{classes}
\keyword{models}