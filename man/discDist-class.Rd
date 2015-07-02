\name{discDist-class}
\docType{class}
\alias{discDist-class}
%\alias{sampleObs,discDist,numeric-method}

\title{Class "discDist"}
\description{Class for discrete probability distributions.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("discDist", alpha, prob, dstr)}.
Either \code{alpha} and \code{prob} or \code{dstr} have to be provided.
The resulting object represents a discrete probability distribution over alphabet \code{alpha}
where the probabilities for individual symbols are given by \code{prob}.
}
\section{Slots}{
	 \describe{
    \item{\code{alpha}:}{A character vector containing all symbols of the alphabet.}
    \item{\code{prob}:}{Numeric vector of same length as \code{alpha} with probabilities
    for each symbol.}
    \item{\code{dstr}:}{A \code{list} of the form \code{symbol = probability}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{dist}"}, directly.
}
\section{Methods}{
  \describe{
    \item{sampleObs}{\code{signature(dist = "discDist", size = "numeric")}: sample observations from \code{alpha}
    with probability \code{prob}. }
	 }
}

\author{Peter Humburg}

\seealso{
	\code{\linkS4class{contDist}}, \code{\linkS4class{tDist}}
}
\examples{
showClass("discDist")
}
\keyword{classes}
\keyword{distribution}
