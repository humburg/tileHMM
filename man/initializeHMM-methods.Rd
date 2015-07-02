\name{initializeHMM-methods}
\docType{methods}
\alias{initialize,contHMM-method}
\alias{new,contHMM-method}

\title{Generate Objects of Class 'hmm'}
\description{
 Objects of class \code{\linkS4class{contHMM}} are created by calling 
 \code{new("contHMM", transition, emission, init)}.
}
\section{Methods}{
\describe{

\item{.Object = "contHMM"}{Creates an object of class \code{contHMM}}
}
}

\usage{
\S4method{new}{contHMM}(Class, transition, emission, init = NULL)
}
\arguments{
	\item{Class}{A character string specifying the class of the object to be created. In this case "contHMM".}
	\item{transition}{A list with one object of class \code{discDist} for each state of the model. These are used as transition probabilities.}
	\item{emission}{A list with one object of class \code{contDist} for each state of the model. These are used as emission distributions.}
	\item{init}{An object of class \code{discDist} or \code{NULL}. This is used as the initial state distribution.}
}
\details{
	The length of all arguments (except \code{Class}) has to be equal. If \code{init} is \code{NULL} a uniform initial 
	state distribution is used.
}

\seealso{
	\code{\linkS4class{contHMM}}, \code{\linkS4class{discDist}}, \code{\linkS4class{contDist}}
}
\author{Peter Humburg}
\keyword{methods}
