\name{generate.data}
\alias{generate.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate Simulated Dataset}
\description{
  Generate simulated data based on real data and the results of a previous analysis.
}
\usage{
generate.data(data, group, pos.range = c(1, 10), 
    num.seq = 100, gap = 35, split.gap = 1000, min.len = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{data.frame} with information about genomic coordinates of probes 
  (chromosome and position) in the first two columns. Subsequent columns contain probe 
  measurements of individual samples.}
  \item{group}{Information that can be used to assign probes to one of two classes. Either a 
  logical vector or the name of a GFF file. In the later case all probes in annotated regions are 
  considered to be \sQuote{positive}.}
  \item{pos.range}{Indicates how many positive regions should be generated for each 
  observation sequence. The actual number for each sequence is sampled uniformly from 
  the indicated range of values.}
  \item{num.seq}{Number of observation sequences to generate.}
  \item{gap}{Gap between probes. Used to generate artificial probe coordinates.}
  \item{split.gap}{Gap between sequences.}
  \item{min.len}{Minimum number of probes per region.}
}

\value{
  A list with components
  \item{observation}{A \code{data.frame} with the same format as \code{data}.}
  \item{regions}{A \code{list} of state sequences.}
}

\author{Peter Humburg}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{datagen}
