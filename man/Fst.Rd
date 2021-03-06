\name{Fst}
\alias{Fst}
\title{
Calculate fixation indices
}
\description{
This function calculates the fixation index Fst for each SNP, together
with its weight in the overall estimate (as used by the Internation
HapMap Consortium). 
}
\usage{
Fst(snps, group, pairwise=FALSE)
}
\arguments{
  \item{snps}{an object of class \code{SnpMatrix} or
    \code{XSnpMatrix} containing the  SNP data}
  \item{group}{a factor (or object than can be coerced into a factor),
    of length equal to the number of rows of \code{snps}, giving the
    grouping or rows for which the Fst is to be calculated}	 
  \item{pairwise}{if \code{TRUE}, the within-group variances are weighted 
    according to the number of possible within-group pairwise comparisons 
    of chromosomes. If \code{FALSE}, the default value, 
    weights are simply the number of chromosomes in each group.}
}
\details{
See vignette.
}
\value{
  A list:
  \describe{
    \item{Fst}{Fst values for each SNP}
    \item{weight}{The weights for combining these into a single index}
    }
}
\note{Uncertain genotypes are treated as missing
}
\author{David Clayton \email{dc208@cam.ac.uk}
}
\examples{

## Analysis of some HapMap data

data(for.exercise)
f <- Fst(snps.10, subject.support$stratum)
weighted.mean(f$Fst, f$weight)

}
\keyword{univar}

