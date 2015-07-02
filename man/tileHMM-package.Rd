\name{tileHMM-package}
\alias{tileHMM-package}
\alias{tileHMM}
\docType{package}
\title{
Hidden Markov Models for ChIP-on-Chip Analysis
}
\description{
This package provides methods and classes to build HMMs that are suitable for the analysis
of ChIP-on-chip data. The provided parameter estimation methods include the Baum-Welch algorithm 
and Viterbi training as well as a combination of both. The latter provides results identical to
the Baum-Welch algorithm but is conciderably faster.
}
\details{
\tabular{ll}{
Package: \tab tileHMM\cr
Type: \tab Package\cr
Version: \tab 1.0-5\cr
Date: \tab 2012-03-19\cr
License: \tab GPL (>= 2)\cr
}
Hidden Markov models are represented as objects of class \code{\linkS4class{hmm}} or derived classes. Function
\code{\link{getHMM}} provides an easy to use interface to create \code{\linkS4class{contHMM}} objects with emission
distributions of class \code{\linkS4class{tDist}} from a set of parameters. Function \code{\link{hmm.setup}} can
be used to create HMMs with initial parameter estimates obtained from data.

To optimise initial parameter estimates \link[=viterbiTraining]{Viterbi training} and the 
\link[=baumWelch]{Baum-Welch} algorithm are provided by this package. Function \code{\link{viterbiEM}} provides
a convenient way to use a combination of both methods.
}
\author{
Peter Humburg

Maintainer: Peter Humburg \email{Peter.Humburg@well.ox.ac.uk}
}
\references{
Humburg, P. and Bulger, D. and Stone, G.
Parameter estimation for robust HMM analysis of ChIP-chip data. BMC Bioinformatics 2008, 9:343.
}

\keyword{ package }
\keyword{models}
\keyword{optimize}

\seealso{
Classes provided by this package:
\code{\linkS4class{hmm}}, \code{\linkS4class{contHMM}}, \code{\linkS4class{dist}}, \code{\linkS4class{discDist}}, 
\code{\linkS4class{contDist}}, \code{\linkS4class{tDist}}

Important functions:
\code{\link{viterbiEM}}, \code{\link{baumWelch}}, \code{\link{viterbiTraining}}
}
