\name{families}
\alias{families}
\alias{genotypes}
\alias{pedData}
\docType{data}
\title{Test data for family association tests}
\description{
  These data started life as real data derived from an affected sibling
  pair study of type 1 diabetes. However, original subject and SNP
  identidiers have been replaced by randomly chosen ones.
}
\usage{data(families)}
\format{
  There are two objects in the loaded data file:
  \itemize{
    \item \code{genotypes}: An object of class
    \code{"\link[=SnpMatrix-class]{snp.matrix}"} containing the SNP
    genotype data for both parents and affected offspring
    \item \code{pedData}: A data frame containing the standard six fields
    for a \emph{LINKAGE} pedfile. The are named \code{familyid},
    \code{member}, \code{father}, \code{mother} \code{sex}, and
    \code{affected}
  }
  The two objects are linked by common row names.
}
\details{
Coding in the \code{pedData} frame is as in the \emph{LINKAGE} package,
except that missing data are coded \code{NA} rather than zero
}
\examples{
data(families)
summary(genotypes)
summary(pedData)
}
\keyword{datasets}
