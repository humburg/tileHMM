\name{getHMM}
\alias{getHMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create HMM from Parameter Values}
\description{
  Create a two state HMM with t emission distributions from a list of parameters.
}
\usage{
getHMM(params, snames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{ A list with components mu, sigma and nu, each a vector with two elements.
   They provide values for the location parameter, scale parameter and degrees of
   freedom for both states. Component a is a vector of length two providing the
   off diagonal elements of the transition probability matrix.}
  \item{snames}{Two character strings that should be used as state names.}
}
\value{
  Object of class contHMM with the parameters provided.
}

\author{Peter Humburg}

\seealso{\code{\link{hmm.setup}}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.035, 0.01)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
hmm <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df), 
    state.names)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
