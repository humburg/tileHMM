\name{tDist-class}
\docType{class}
\alias{tDist-class}

\title{Class "tDist"}
\description{This class provides objects representing scaled \emph{t} distributions.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tDist", mean, var, df)}.
	\describe{
      \item{\code{mean}:}{Location parameter.}
      \item{\code{var}:}{Scale parameter.}
      \item{\code{df}:}{Degrees of freedom.}
	}
	The distribution has density
	\deqn{%
		f(x; \texttt{mean}, \texttt{var}, \texttt{df}) = \frac{\Gamma\left(\frac{\small\texttt{df}+1}{2}\right)}{%
		\sqrt{\small\texttt{var}\times\pi\times\small\texttt{df}}\,\Gamma\left(\frac{\small\texttt{df}}{2}\right)} \times %
		\frac{1}{\sqrt{\left(1 + \frac{(x-\small\texttt{mean})^2}{\small\texttt{var}\times\small\texttt{df}}\right)^{\small\texttt{df}+1}}}}{%
		f(x; mean, var, df) = gamma((df+1)/2)/\{sqrt[var * pi * df (1 + (x - mean)^2 / (var * df))^(df + 1)] * gamma(df/2)\}}
}
\section{Slots}{
	 \describe{
    \item{\code{components}:}{Object of class \code{"matrix"} with columns \code{weight}, \code{mean}, 
    \code{variance}, and \code{df}.  For objects of class \code{"tDist"} this matrix has only one row 
    and \code{weight} is always 1.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{contDist}"}, directly.
Class \code{"\linkS4class{dist}"}, by class "contDist", distance 2.
}
\section{Methods}{
  \describe{
    \item{sampleObs}{\code{signature(dist = "tDist", size = "numeric")}: Sample \code{size} observations from \code{dist}.}
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Peter Humburg}

\seealso{
	\code{\linkS4class{contHMM}}
}
\examples{
showClass("tDist")
}
\keyword{classes}
\keyword{distribution}
