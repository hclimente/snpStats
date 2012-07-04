\name{ibsDist}
\alias{ibsDist}
\title{Distance matrix based on identity by state (IBS)}
\description{
Expresses a matrix of IBS counts (see \code{\link{ibsCount}}) as a
distance matrix. The distance between two samples is returned as the
proportion of allele comparisons which are \emph{not} IBS.
}
\usage{
ibsDist(counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{A matrix of IBS counts as produced by the function
    \code{\link{ibsCount}}}
}
\value{
  An object of class \code{"dist"} (see \code{\link[stats]{dist}})
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{ibsCount}}, \code{\link[stats]{dist}}}
\examples{
data(testdata)
ibs <- ibsCount(Xchromosome)
distance <- ibsDist(ibs)
}
\keyword{cluster}

