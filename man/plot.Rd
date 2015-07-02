\name{plot}
\alias{plot.contDist}
\alias{plot.contHMM}
\alias{plot.tDist}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting of "contDist" Objects}
\description{
  Functions for plotting of contDist objects, either directly or as part of contHMM objects.
}
\usage{
\method{plot}{contDist}(x, step.size = 0.01, new.plot = TRUE, weight = 1, ...)
\method{plot}{tDist}(x, step.size = 0.01, new.plot = TRUE, weight = 1, ...)
\method{plot}{contHMM}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class contDist or contHMM.}
  \item{step.size}{Numeric value indicating the distance between points at which the density function
    is evaluated.}
  \item{new.plot}{If this is TRUE (the default) a new plot is created, otherwise graph of the density
   function is added to the current plot.}
  \item{weight}{Weighting factor. The density function will be scaled by this factor.
   This is useful when plotting mixture components. }
  \item{\dots}{Additional arguments to be passed to plot}
}

\value{
These functions are called for their side effect.
}
\author{Peter Humburg}

\seealso{\linkS4class{contDist}, \linkS4class{contHMM}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.1, 0.02)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
model <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df),
    state.names)
    
## plot emission distributions of HMM
plot(model)

## now plot density functions separately
par(mfrow=c(1,1))
plot(model@emission$one)
plot(model@emission$two, new.plot=FALSE, lty=2)
legend("topleft", legend=states(model), lty=1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
