\name{shrinkt.st}
\alias{shrinkt.st}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate 'Shrinkage t' Statistic}
\description{
Calculate the \sQuote{shrinkage \emph{t}} statistic (Opgen-Rhein and Strimmer, 2007).
}
\usage{
shrinkt.st(X, L, h0.mean = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Data matrix. Rows correspond to variables (probes) and columns to samples.}
  \item{L}{Design vector. Indicating treatment (1) and control (2) samples. If no control
    samples are given a one sample \emph{t} test is carried out.}
  \item{h0.mean}{If no control samples are provided the treatment mean is compared to this
    value in a one sample \emph{t} test.}
  \item{\dots}{Further arguments to be passed to and from other methods.}
}

\value{
 An object of class "\code{shrinkage}" containing the test statistics.
}
\references{
Opgen-Rhein, R., and K. Strimmer 2007. 
Accurate ranking of differentially expressed genes by a distribution-free shrinkage 
approach.   \emph{Statist. Appl. Genet. Mol. Biol.} \strong{6}:9. 
\url{http://www.bepress.com/sagmb/vol6/iss1/art9/}
}
\author{Peter Humburg}
\note{
This method uses the \code{\link[st]{shrinkt.stat}} method from package \code{st}.
}
\seealso{\code{\link[st]{shrinkt.stat}}}
\examples{
## generate some data
X <- matrix(nrow = 100, ncol = 6)
X[ , 1:3] <- apply(X[ , 1:3], 1, rnorm, 
    mean = rnorm(1, mean = 2, sd = 2), sd = rchisq(1, df = 2))
X[ , 4:6] <- apply(X[ , 4:6], 1, rnorm, 
    mean = rnorm(1, mean = 0, sd = 1), sd = rchisq(1, df = 2))
L <- c(1, 1, 1, 2, 2, 2)

## compute shrinkage t statistic
st.stat <- shrinkt.st(X, L)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{ utilities }
