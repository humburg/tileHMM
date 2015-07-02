\name{contDist-class}
\docType{class}
\alias{contDist-class}
\alias{show,contDist-method}
%\alias{.baumWelchEmission,contHMM,contDist,list-method}

\title{Class "contDist"}
\description{Class for continuous probability distributions.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("contDist", weight, center, disp)}.
These objects represent a continuous probability distribution as a mixture of one or more
elementary distributions.
}
\section{Slots}{
	 \describe{
    \item{\code{components}:}{Object of class \code{"matrix"} with one row for each mixture component
    and a column for each of mixture weight, location parameter and dispersion parameter.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{dist}"}, directly.
}
\section{Methods}{
\describe{
\item{show}{\code{signature(object = "contDist")}: Print summary of contDist object.}
}
}

\author{Peter Humburg}

\seealso{
	\code{\linkS4class{tDist}}, \code{\linkS4class{discDist}}
}
\examples{
showClass("contDist")
}
\keyword{classes}
\keyword{distribution}
