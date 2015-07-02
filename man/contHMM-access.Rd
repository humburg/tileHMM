\name{contHMM-access}
\alias{contHMM-access}
\alias{[.contHMM}
\alias{length,hmm-method}
\title{ Accessing Objects of Class "contHMM"}
\description{
 Access to model parameters and densities of emission distributions.
}
\usage{
\method{[}{contHMM}(x, i, j, transition = TRUE, log = FALSE, sum = TRUE, ...)
\S4method{length}{hmm}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{contHMM}}
  \item{i}{State for which parameter values should be retrieved. This can either a
     numeric value giving the state index or a character string with the state name.}
  \item{j}{Second index identifying parameter (see Details).}
  \item{transition}{Logical indicating whether transition probabilities or density values should be returned.}
  \item{log}{Logical indicating whether values should be log transformed before they are returned.}
  \item{sum}{Logical indicating whether densities of mixture components should be summed up.
     This is ignored if \code{transition = TRUE}.}
  \item{\dots}{Futher arguments to be passed to and from other methods.}
}
\details{
 The \sQuote{[} function allows access to the transition probability matrix of the model as well as the emission 
 distributions. If \code{transition = TRUE} the transition probability matrix is accessed. In this case
 \code{i} and \code{j} identify rows and columns of the matrix respectively. Both can be given as either 
 numeric index or name of the respective states. Either or both of \code{i} and \code{j} may be missing
 to indicate that an entire row or column should be selected.
 
 If \code{transition = FALSE} the emission distribution of state \code{i} is accessed instead. In this case
 the density function is evaluated at point \code{j}. 
}
\value{
  For \sQuote{[} either a subset of the transition probability matrix of \code{x} or the probability density of state \code{i}
  evaluated at point \code{j} (see Details).
  For \sQuote{length} the number of states in the model.
}
\author{Peter Humburg}

\seealso{\code{\linkS4class{contHMM}}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.1, 0.02)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
hmm <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df), 
    state.names)

## number of states in the model
length(hmm)

## transition probability from state 'one' to state 'two'
hmm["one", "two"]
## or equivalently
hmm[1, 2]

## get the transition probability matrix
hmm[ , ]

## evaluate emission distribution function of state 'one' at 0
hmm["one", 0, transition = FALSE]

## again, this time using log transformation
hmm["one", 0, transition = FALSE, log = TRUE]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
