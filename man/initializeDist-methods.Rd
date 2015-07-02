\name{initializeDist-methods}
\docType{methods}
\alias{initialize,contDist-method}
\alias{initialize,discDist-method}
\alias{initialize,tDist-method}
\alias{new,contDist-method}
\alias{new,discDist-method}
\alias{new,tDist-method}


\title{Generating Objects of Class 'dist'}
\description{
 Objects representing probability distributions for the use with objects of class
 \code{hmm} are created by calls to \code{new(Class, ...)}, where \code{Class} is
 \code{"contDist"}, \code{"tDist"} or \code{"discDist"}. 
}
\section{Methods}{
\describe{

\item{.Object = "contDist"}{Creates an object of class \code{contDist}, 
representing a continuous probability distribution. This is used as a super class 
for specific distributions and should not be used directly.}

\item{.Object = "tDist"}{Creates an object of class \code{tDist}, 
representing a t distribution with location and scale parameter.}

\item{.Object = "discDist"}{Creates an object of class \code{discDist}, 
representing a discrete probability distribution.}
}
}

\usage{
\S4method{new}{tDist}(Class, mean=0, var=1, df=3)
\S4method{new}{discDist}(Class, alpha=character(0), prob=numeric(0), dstr=list())
}
\arguments{
	\item{Class}{Character string specifying the class of the object to be created.}
	\item{mean}{The location parameter of the \emph{t} distribution.}
	\item{var}{The scale parameter of the \emph{t} distribution.}
	\item{df}{The degrees of freedom of the \emph{t} distribution.}
	\item{alpha}{A character vector with names for all events.}
	\item{prob}{A numeric vector providing the probabilities of events.}
	\item{dstr}{A list with each entry providing the probability of one event.}
}

\details{
	To generate objects of class \code{discDist} either \code{alpha} and \code{prob} 
	or \code{dstr} have to be provided. If \code{alpha} and \code{prob} are provided 
	they have to be the same length.
}

\seealso{\code{\linkS4class{discDist}}, \code{\linkS4class{contDist}}, \code{\linkS4class{tDist}}}
\author{Peter Humburg}

\keyword{methods}
