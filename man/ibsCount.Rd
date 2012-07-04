\name{ibsCount}
\alias{ibsCount}
\title{Count alleles identical by state}
\description{
  This function counts, for all pairs of subjects and across all SNPs,
  the total number of alleles which are identical by state (IBS)
}
\usage{
ibsCount(snps, uncertain = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{An input object of class \code{"SnpMatrix"} or
    \code{"XSnpMatrix"}}
  \item{uncertain}{If \code{FALSE}, uncertain genotypes are
    ignored. Otherwise contributions are weighted by posterior probabilities}
}
\details{
For each pair of subjects the function counts the total number of
alleles which are IBS. For autosomal SNPs, each locus
contributes 4 comparisons, since each subject carries two
copies. For SNPs on the X chromosome, the number of comparisons is also
4 for female:female comparisons, but is 2 for female:male 
and 1 for male:male comparisons.
}
\value{
  If there are \var{N} rows in the input matrix, the function returns an
  \var{N}*\var{N} matrix. The upper triangle contains the total number
  of comparisons and the lower triangle contains the number of these
  which  are IBS. The diagonal contains the number of valid calls for
  each subject. 
}
\note{
In genome-wide studies, the SNP data will usually be held as a series of
objects (of
class \code{"SnpMatrix"} or \code{"XSnpMatrix"}), one 
per chromosome. Note that the matrices
produced by applying the \code{ibsCount} function to each object in turn
can be added to yield the genome-wide result.
}
\author{David Clayton \email{dc208@cam.ac.uk} }
\seealso{\code{\link{ibsDist}} which calculates a distance matrix based
  on proportion of alleles which are IBS}
\examples{
data(testdata)
% taking a subset of Autosomes - it's a bit big for testing basic functionality
ibs.A <- ibsCount(Autosomes[,1:100])
ibs.X <- ibsCount(Xchromosome)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

