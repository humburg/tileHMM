\name{remove.short}
\alias{remove.short}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-Processing of "tileHMM" Results}
\description{
  Remove short regions that are likely to be spurious.
}
\usage{
remove.short(regions, post, probe.pos, min.length = 1000, 
    min.score = 0.8, summary.fun = mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regions}{A matrix with information about the location of enriched regions.}
  \item{post}{A numeric vector with the posterior probability of ChIP enrichment for each probe.}
  \item{probe.pos}{A data frame with columns \sQuote{chromosome} and \sQuote{position} providing 
      genomic coordinates for each probe.}
  \item{min.length}{Minimum length of enriched regions (see Details).}
  \item{min.score}{Minimum score for enriched regions (see Details).}
  \item{summary.fun}{Function used to summarise posterior probe probabilities into region scores.}
}
\details{
  All regions that are shorter than \code{min.length} \strong{and} have a score of less than 
  \code{min.score} will be removed. To filter regions based on only one of these values set the 
  other one to 0.
  
  Region scores are calculated based on posterior probe probabilities. The summary function used
  should accept a single numeric argument and return a numeric vector of length 1. If the
  probabilities in \code{post} are log transformed they will be transformed back to linear
  space before they are summarised for each region.
}
\value{
  A matrix with two rows and one column for each remaining region.
}

\author{Peter Humburg}

\seealso{\code{\link{region.position}}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.1, 0.1)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
model <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df), 
    state.names)

## obtain observation sequence from model
obs <- sampleSeq(model, 500)

## make up some genomic probe coordinates
pos <- data.frame(chromosome = rep("chr1", times = 500), 
    position = seq(1, 18000, length = 500))

## calculate posterior probability for state "one"
post <- posterior(obs, model, log=FALSE)

## get sequence of individually most likely states
state.seq <- apply(post, 2, which.max)
state.seq <- states(model)[state.seq]

## find regions attributed to state "one"
reg.pos <- region.position(state.seq, region="one")

## remove short and unlikely regions
reg.pos2 <- remove.short(reg.pos, post, pos, min.length = 200, 
    min.score = 0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

