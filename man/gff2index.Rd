\name{gff2index}
\alias{gff2index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract Probe Calls from GFF File}
\description{
Creates a logical vector indicating probes in annotated regions from information
in a GFF formatted file.
}
\usage{
gff2index(gff, pos)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gff}{GFF file (see Details).}
  \item{pos}{A \code{data.frame} with chromosome names in its first column and 
  probe positions in the second column.}
}
\details{
  The GFF file can be provided in several ways. Either as \code{data.frame}, as 
  a character string providing the name of a GFF file or as a \code{connection} object 
  pointing to a GFF file.
}
\value{
  A logical vector with one entry for each probe in \code{pos}. \code{TRUE} indicates probes that
  are inside a region that is annotated in the provided GFF file. 
}
\references{GFF specification: \url{http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml}}
\author{Peter Humburg}

\seealso{\code{\link{region.position}} and \code{\link{reg2gff}} for the reverse operation.}
\examples{
## create two state HMM with t distributions
state.names <- c("one","two")
transition <- c(0.1, 0.02)
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

## find regions attributed to state "two"
reg.pos <- region.position(vit.res$stateSeq, region="two")

## calculate posterior probability for state "two"
post <- posterior(obs, model, log=FALSE)[1,]

## create gff annotations
gff <- reg2gff(reg.pos, post, pos)

## extract probe calls from gff annotation
probe.idx <- gff2index(gff, pos)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
