\name{xxt}
\alias{xxt}
\title{X.X-transpose for a standardized SnpMatrix}
\description{
The input SnpMatrix is first standardized by subtracting the mean (or
stratum mean) from
each call  and
dividing by the expected standard deviation under Hardy-Weinberg
equilibrium.
It is then post-multiplied by its transpose. This is a
preliminary step in the computation of principal components.
}
\usage{
xxt(snps, strata = NULL, correct.for.missing = FALSE, lower.only = FALSE,
  uncertain = FALSE)
}
\arguments{
  \item{snps}{The input matrix, of type \code{"SnpMatrix"}}
  \item{strata}{A \code{factor} (or an object which can be coerced into
    a  \code{factor}) with length equal to the number of rows of
    \code{snps} defining stratum membership}
  \item{correct.for.missing}{If \code{TRUE}, an attempt is made to
    correct for the effect of missing data by use of inverse probability
    weights. Otherwise, missing observations are scored zero in the
    standardized matrix}
  \item{lower.only}{If \code{TRUE}, only the lower triangle of the
    result is returned and the upper triangle is filled with
    zeros. Otherwise, the complete symmetric matrix is returned}
  \item{uncertain}{If \code{TRUE}, uncertain genotypes are replaced by
    posterior expectations. Otherwise these are treated as missing
    values}
}
\details{
  This computation forms the first step of the calculation of principal
  components for genome-wide SNP data. As pointed out by Price et al.
  (2006), when the data matrix has more rows than columns it is
  most efficient to calculate the eigenvectors of
  \var{X}.\var{X}-transpose, where \var{X} is a 
  \code{SnpMatrix} whose columns have been  standardized to zero mean and
  unit variance. For autosomes, the genotypes are given codes 0, 1 or 2
  after subtraction of the mean, 2\var{p}, are divided by the standard
  deviation 
  sqrt(2\var{p}(1-\var{p})) (\var{p} is the estimated allele
  frequency). For SNPs on the X chromosome in male subjects,
  genotypes are coded 0 or 2. Then 
  the mean is still 2\var{p}, but the standard deviation is 
  2sqrt(\var{p}(1-\var{p})). If the \code{strata} is supplied, a
  stratum-specific estimate value for \var{p} is used for
  standardization. 
  
  Missing observations present some difficulty. Price et al. (2006)
  recommended replacing missing observations by their means, this being
  equivalent to replacement by zeros in the standardized matrix. However
  this results in a biased estimate of the complete data
  result. Optionally this bias can be corrected by inverse probability
  weighting. We assume that the probability that any one call is missing
  is small, and can be predicted by a multiplicative model with row
  (subject) and column (locus) effects. The estimated probability of a
  missing value in a given row and column is then given by
  \eqn{m = RC/T}, where \var{R} is the row total number of
  no-calls, \var{C} is the column total of no-calls, and \var{T} is the
  overall total number of no-calls. Non-missing contributions to
  \var{X}.\var{X}-transpose are then weighted by \eqn{w=1/(1-m)} for
  contributions to the diagonal elements, and products of the relevant
  pairs of weights for contributions to off--diagonal elements.
}
\value{
A square matrix containing either the complete X.X-transpose matrix, or
just its lower triangle
}
\references{
Price et al. (2006) Principal components analysis corrects for
stratification in genome-wide association studies. \emph{Nature Genetics},
\bold{38}:904-9
}
\note{
In genome-wide studies, the SNP data will usually be held as a series of
objects (of
class \code{"SnpMatrix"} or\code{"XSnpMatrix"}), one per chromosome.
Note that the  \var{X}.\var{X}-transpose matrices
produced by applying the \code{xxt} function to each object in turn
can be added to yield the genome-wide result.

If the matrix is converted to a correlation matrix by pre- and
post-multiplying by the sqrt of the inverse of its diagonal, then this
is an unbiased estimate of twice the kinship matrix.
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\section{Warning}{
  The correction for missing observations can result in an output matrix
  which is not positive semi-definite. This should not matter in the
  application for which it is intended
  }
\examples{
# make a SnpMatrix with a small number of rows
data(testdata)
small <- Autosomes[1:100,]
# Calculate the X.X-transpose matrix
xx <- xxt(small, correct.for.missing=TRUE)
# Calculate the principal components
pc <- eigen(xx, symmetric=TRUE)$vectors
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}
\keyword{multivariate}

