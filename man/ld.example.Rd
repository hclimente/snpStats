\name{ld.example}
\docType{data}
\alias{ld.example}
\alias{ceph.1mb}
\alias{yri.1mb}
\alias{support.ld}
\title{Datasets to illustrate calculation of linkage disequilibrium statistics}
\description{
  This R data file contains data from the International HapMap project,
  concerning 603 SNPs spanning a one megabase region on chromosome 22,
  in a sample of Europeans and a sample of Africans
}
\format{
  There are three objects in the file:
  \itemize{
    \item \code{ceph.1m}: A \code{snpMatrix} object containing the
    European genotype data
    \item \code{yri.1m}: A \code{snpMatrix} object containing the
    African genotype data
    \item \code{support.ld}: A dataframe containing details (chromosome
    position etc.} of the 603 SNPs
}
\source{
  \url{http://hapmap.ncbi.nlm.nih.gov}
}
\references{
  The International HapMap Consortium. The International HapMap Project.
  \emph{Nature} \bold{426}:789-796 (2003)
}
\keyword{datasets}
