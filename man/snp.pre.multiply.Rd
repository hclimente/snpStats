\name{snp.pre.multiply}
\alias{snp.pre.multiply}
\alias{snp.post.multiply}
\title{Pre- or post-multiply a SnpMatrix object by a general matrix}
\description{
  These functions first standardize the input \code{SnpMatrix} in the
  same way as does the function \code{\link{xxt}}. The standardized
  matrix is then either pre-multiplied (\code{snp.pre.multiply}) or
  post-multiplied (\code{snp.post.multiply}) by a general matrix. Allele
  frequencies for standardizing the input SnpMatrix may be supplied
  but, otherwise, are calculated from the input SnpMatrix
}
\usage{
snp.pre.multiply(snps, mat, frequency=NULL, uncertain = FALSE)
snp.post.multiply(snps, mat, frequency=NULL, uncertain = FALSE)
}
\arguments{
  \item{snps}{An object of class \code{"SnpMatrix"} or \code{"XSnpMatrix"}}
  \item{mat}{A general (numeric) matrix}
  \item{frequency}{A numeric vector giving the allele (relative)
    frequencies to be used for standardizing the columns of \code{snps}.
    If \code{NULL}, allele frequencies will be calculated
    internally. Frequencies should refer to the second (\code{B}) allele
  }
  \item{uncertain}{If \code{TRUE}, uncertain genotypes are replaced by
    posterior expectations. Otherwise these are treated as missing values}
}
\details{
The two matrices must be conformant, as with standard matrix
multiplication. The main use envisaged for these functions is the
calculation of factor loadings in principal component analyses of large
scale SNP data, and the application of these loadings to other
datasets. The use of externally supplied allele frequencies for
standardizing the input SnpMatrix is required when applying loadings
calculated from one dataset to a different dataset
}
\value{
  The resulting matrix product
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{xxt}}}
\examples{
##--
##-- Calculate first two principal components and their loading, and verify
##--	
# Make a SnpMatrix with a small number of rows
data(testdata)
small <- Autosomes[1:20,]
# Calculate the X.X-transpose matrix
xx <- xxt(small, correct.for.missing=FALSE)
# Calculate the first two principal components and corresponding eigenvalues
eigvv <- eigen(xx, symmetric=TRUE)
pc <- eigvv$vectors[,1:2]
ev <- eigvv$values[1:2]
# Calculate loadings for first two principal components
Dinv <- diag(1/sqrt(ev))
loadings <- snp.pre.multiply(small,  Dinv \%*\% t(pc))
# Now apply loadings back to recalculate the principal components
pc.again <- snp.post.multiply(small, t(loadings) \%*\% Dinv)
print(cbind(pc, pc.again))
}
\keyword{array}
\keyword{multivariate}

