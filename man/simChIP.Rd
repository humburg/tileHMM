\name{simChIP}
\alias{simChIP}
\docType{data}
\title{Simulated ChIP-on-Chip Data}
\description{
A simulated ChIP-on-chip dataset. This is part of the data used in a simulation study
by Humburg \emph{et al.} (2007), it is based on data published by Zhang \emph{et al.} (2007).
}
\usage{data(simChIP)}
\format{
  A data frame with 210828 probe measurements under three different conditions with four replicates each.
  \describe{
    \item{\code{chromosome}}{The chromosome targeted by probe. Here always \sQuote{chr1}.}
    \item{\code{position}}{Position of central base.}
    \item{\code{input_1}}{Genomic DNA sample 1.}
    \item{\code{input_2}}{Genomic DNA sample 2.}
    \item{\code{input_3}}{Genomic DNA sample 3.}
    \item{\code{input_4}}{Genomic DNA sample 4.}
    \item{\code{h3_1}}{Histone H3 ChIP sample 1.}
    \item{\code{h3_2}}{Histone H3 ChIP sample 2.}
    \item{\code{h3_3}}{Histone H3 ChIP sample 3.}
    \item{\code{h3_4}}{Histone H3 ChIP sample 4.}
    \item{\code{h3k27_1}}{Histone H3K27me3 ChIP sample 1.}
    \item{\code{h3k27_2}}{Histone H3K27me3 ChIP sample 2.}
    \item{\code{h3k27_3}}{Histone H3K27me3 ChIP sample 3.}
    \item{\code{h3k27_4}}{Histone H3K27me3 ChIP sample 4.}
  }
}
\source{
  Humburg, P. and Bulger, D. and Stone, G. 2008
Parameter estimation for robust HMM analysis of ChIP-chip data. \emph{unpublished}
}
\references{
  Zhang, X. and Clarenz, O. and Cokus, S. and Bernatavichute, Y. V. and Goodrich, J. and Jacobsen S. E. 2007
  Whole-genome analysis of histone H3 lysine 27 trimethylation in Arabidopsis. 
  \emph{PLoS Biology}, \bold{5}(5).
}
\examples{
data(simChIP)
\dontrun{
## scatter plots
library(geneplotter)
simChIP[ , 3:14] <- logb(simChIP[ , 3:14], 2)
par(mfrow=c(2,2))
smoothScatter(simChIP[ , 7], simChIP[ , 11], xlab = "H3 sample 1", 
    ylab = "H3K27me3 sample 1")
smoothScatter(simChIP[ , 8], simChIP[ , 12], xlab = "H3 sample 2", 
    ylab = "H3K27me3 sample 2")
smoothScatter(simChIP[ , 9], simChIP[ , 13], xlab = "H3 sample 3", 
    ylab = "H3K27me3 sample 3")
smoothScatter(simChIP[ , 10], simChIP[ , 14], xlab = "H3 sample 4", 
    ylab = "H3K27me3 sample 4")
}
}
\keyword{datasets}
