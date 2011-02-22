\name{single.snp.tests}
\alias{single.snp.tests}
\title{1-df and 2-df tests for genetic associations with SNPs (or
  imputed SNPs)}
\description{
  This function carries out tests for association between phenotype and
  a series of single nucleotide polymorphisms (SNPs), within strata
  defined by a possibly confounding factor. SNPs are
  considered one at a time and both 1-df and 2-df tests are
  calculated. For a binary phenotype, the
  1-df test is the Cochran-Armitage test (or, when
  stratified, the Mantel-extension test). The function will also calculate
  the same tests for SNPs imputed by regression analysis.
}
\usage{
single.snp.tests(phenotype, stratum, data = sys.parent(), snp.data,
   rules=NULL, subset, snp.subset, uncertain = FALSE, score=FALSE)
}
\arguments{
  \item{phenotype}{A vector containing the values of the phenotype}
  \item{stratum}{Optionally, a factor defining strata for the analysis}
  \item{data}{A dataframe containing the \code{phenotype} and
    \code{stratum} data. The row names of this are linked with the row
    names of the \code{snps} argument to establish correspondence of
    phenotype and genotype data. If this argument is not supplied,
    \code{phenotype} and \code{stratum} are evaluated in the calling
    environment and should be in the same order as rows of \code{snps}}
  \item{snp.data}{An object of class \code{"SnpMatrix"} containing the SNP
    genotypes to be tested}
  \item{rules}{An object of class
    \code{"ImputationRules"}. If
    supplied, the rules coded in this object are used, together with
    \code{snp.data}, to calculate tests for imputed SNPs}
  \item{subset}{A vector or expression describing the subset of subjects
    to be used in the analysis. This is evaluated in the same
    environment as the \code{phenotype} and \code{stratum} arguments}
  \item{snp.subset}{A vector describing the subset of SNPs to be
    considered. Default action is to test all SNPs in \code{snp.data}
    or, in imputation mode,  as specified by \code{rules}}
  \item{uncertain}{If TRUE, uncertain genotypes are handled by replacing
    score contributions by their posterior expectations. Otherwise they
    are treated as missing. Setting this option authomatically invokes
    use of \code{robust} variance estimates}
   \item{score}{If \code{TRUE}, the output object will contain, for each
    SNP,  the score vector and its variance-covariance matrix}
}
\details{
  Formally, the test statistics are score tests for generalized linear
  models with canonical link. That is, they are inner products between
  genotype indicators and the deviations of phenotypes from their
  stratum means. Variances (and covariances) are those of the
  permutation distribution obtained by randomly permuting phenotype
  within stratum.

  When the function is used to calculate tests for imputed SNPs, the
  test is still a score test. The score statistics are calculated from
  the expected value, given observed SNPs, of the score statistic if the
  SNP to be tested were itself observed.
  
  The \code{subset} argument can either be a logical vector of length
  equal to the length of the vector of  phenotypes, an integer vector
  specifying positions in the \code{data} frame, or a character vector
  containing names of the selected rows in the \code{data}
  frame. Similarly, the \code{snp.subset} argument can be a logical,
  integer, or character vector.
}
\value{
  An object of class
  \code{"SingleSnpTests"}. 
  If \code{score} is set to \code{TRUE},
  the output object will be of the extended class
  \code{"SingleSnpTestsScore"}
  containing additional slots holding the score statistics and their
  variances (and covariances). This allows meta-analysis using the
  \code{\link{pool}} function.
}
\references{
  Chapman J.M., Cooper J.D., Todd J.A. and Clayton D.G. (2003)
  \emph{Human Heredity}, \bold{56}:18-31.\cr
  Clayton (2008) Testing for association on the X chromosome
  \emph{Biostatistics}, \bold{9}:593-600.)
}
\note{
  The 1 df imputation tests are described by Chapman et al. (2008)
  and the 2 df imputation tests are a simple extension of these.
  The behaviour of this function for objects of class
  \code{XSnpMatrix} is as described by Clayton (2008). Males are
  treated as homozygous females and corrected variance estimates are
  used.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.lhs.tests}}, \code{\link{snp.rhs.tests}},
  \code{\link{impute.snps}}, \code{\link{ImputationRules-class}},
  \code{\link{pool}}, 
	\code{\link{SingleSnpTests-class}},
    \code{\link{SingleSnpTestsScore-class}}
}
\examples{
data(testdata)
results <- single.snp.tests(cc, stratum=region, data=subject.data,
   snp.data=Autosomes, snp.subset=1:10)
print(summary(results))

# writing to an (anonymous and temporary) csv file
csvfile <- tempfile()
write.csv(file=csvfile, as(results, 'data.frame'))
unlink(csvfile)
# QQ plot 
qq.chisq(chi.squared(results, 1), 1)
qq.chisq(chi.squared(results, 2), 2)
}
\keyword{htest}

