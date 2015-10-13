\name{impute.snps}
\alias{impute.snps}
\title{Impute snps}
\description{
Given SNPs stored in an object of class \code{"SnpMatrix"} or
\code{"XSnpMatrix"} and a set of imputation rules in an object of
class \code{"ImputationRules"}, this function calculates imputed values.
}
\usage{
impute.snps(rules, snps, subset = NULL, as.numeric = TRUE)
}
\arguments{
  \item{rules}{The imputation rules; an object of class
    \code{"ImputationRules"}}
  \item{snps}{The object of class \code{"SnpMatrix"} or
    \code{"XSnpMatrix"} containing the observed SNPs}
  \item{subset}{A vector describing the subset of subjects
    to be used. If \code{NULL} (default), then use all subjects}
  \item{as.numeric}{If \code{TRUE}, the output is a numeric matrix
    containing posterior expectations of the imputed SNPs. Otherwise the
    output matrix is of the same class as \code{snps} and contains
    uncertain genotype calls}
}
\value{
  A matrix with imputed SNPs as columns. The imputed values are the
  estimated expected values of each SNP when coded 0, 1 or 2. 
}
\references{
  Wallace, C. et al. (2010) \emph{Nature Genetics}, \bold{42}:68-71
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{snp.imputation}}}
\examples{
# Remove 5 SNPs from a datset and derive imputation rules for them
data(for.exercise)
sel <- c(20, 1000, 2000, 3000, 5000)
to.impute <- snps.10[,sel]
impute.from <- snps.10[,-sel]
pos.to <- snp.support$position[sel]
pos.fr <- snp.support$position[-sel]
imp <- snp.imputation(impute.from, to.impute, pos.fr, pos.to)
# Now calculate the imputed values
imputed <- impute.snps(imp, impute.from)
}
\keyword{models}
\keyword{regression}

