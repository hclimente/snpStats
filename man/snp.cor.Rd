\name{snp.cor}
\alias{snp.cor}
\title{Correlations with columns of a SnpMatrix}
\description{
This function calculates Pearson correlation coefficients  between columns of a
\code{SnpMatrix} and columns of an ordinary matrix. The two matrices
must have the same number of rows. All valid pairs are used in the
computation of each correlation coefficient.
} 
\usage{
snp.cor(x, y, uncertain = FALSE)
}
\arguments{
  \item{x}{An \var{N} by \var{M} \code{SnpMatrix}}
  \item{y}{An \var{N} by \var{P} general matrix}
  \item{uncertain}{If \code{TRUE}, uncertain genotypes are replaced by
    posterior expectations. Otherwise these are treated as missing
    values}
}
\details{
  This can be used together with \code{\link{xxt}} and
  \code{\link[base]{eigen}} to calculate standardized loadings in the principal
  components
}
\value{
An  \var{M} by \var{P} matrix of correlation coefficients
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\note{
  This version cannot handle X chromosomes
}
\seealso{\code{\link{xxt}}}
\examples{
# make a SnpMatrix with a small number of rows
data(testdata)
small <- Autosomes[1:100,]
# Calculate the X.X-transpose matrix
xx <- xxt(small, correct.for.missing=TRUE)
# Calculate the principal components
pc <- eigen(xx, symmetric=TRUE)$vectors
# Calculate the loadings in first 10 components */
loadings <- snp.cor(small, pc[,1:10])
}
\keyword{array}
\keyword{multivariate}


