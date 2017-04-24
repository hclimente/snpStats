\name{random.snps}
\alias{random.snps}
\title{
Generate random SnpMatrix 
}
\description{
This function is purely for testing purposes. It can generate SnpMatrix
objects which contain more than 2^31-1 elements.
}
\usage{
random.snps(nrows, ncols)
}
\arguments{
  \item{nrows}{The number of rows to be generated}
  \item{ncols}{The number of columns to be generated}
}
\details{
All SNPs should be in Hardy-Weinberg equilibrium with an allele
frequency of 0.5.

Note that, although the total number of elements can exceed 2^31-1, the
numbers of rows and columns are still subject to this limit.
}
\examples{
x <- random.snps(100,10)
col.summary(x)
}
\keyword{utilities}
