\name{reg2gff}
\alias{reg2gff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converting Information about Enriched Regions into GFF Format}
\description{
Creates a GFF formatted file with information about enriched regions.
}
\usage{
reg2gff(regions, post, probe.pos, file = NULL, score.fun = mean, 
    source = "tHMM", feature.type = "posterior_prob", 
    class = "ChIP_region", name = "tHMM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regions}{A matrix indicating enriched regions (see \code{\link{region.position}}).}
  \item{post}{A vector with posterior probabilities for each probe.}
  \item{probe.pos}{A \code{data.frame} with columns "\code{chromosome}" and "\code{position}" providing genomic
  coordinates for each probe.}
  \item{file}{Name of GFF file to create.}
  \item{score.fun}{Function used to calculate score of enriched regions.}
  \item{source}{Entry for \sQuote{source} field of GFF file.}
  \item{feature.type}{Entry for \sQuote{feature} field of GFF file.}
  \item{class}{Class of documented feature. This is used in the \sQuote{attribute} field together with \code{name}.}
  \item{name}{Name of documented feature. This is used in the \sQuote{attribute} field together with \code{class}.}
}
\details{
\code{post} should provide scores for each probe. These scores are then summarised for each region by applying
\code{score.fun} to the probe scores in each region. The default for \code{score.fun} is \code{mean} but
any function that accapts a single numeric vector as its argument and returns a scalar can be used.
}
\value{
A GFF formated \code{data.frame} is returned invisibly. Usually this function is called for
its side effect.
}
\references{GFF specification: \url{http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml}}
\author{Peter Humburg}

\seealso{\code{\link{region.position}}}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.1, 0.2)
location <- c(1, 2)
scale <- c(1, 1)
df <- c(4, 6)
model <- getHMM(list(a=transition, mu=location, sigma=scale, nu=df),
    state.names)

## obtain observation sequence from model
obs <- sampleSeq(model, 100)

## make up some genomic probe coordinates
pos <- data.frame(chromosome = rep("chr1", times = 100), 
    position = seq(1, 4000, length = 100))

## compute most likely state sequence for obs
vit.res <- viterbi(model, obs)

## find regions attributed to state "one"
reg.pos <- region.position(vit.res$stateSeq, region="one")

## calculate posterior probability for state "one"
post <- posterior(obs, model, log=FALSE)[1,]

## create gff annotations
gff <- reg2gff(reg.pos, post, pos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}