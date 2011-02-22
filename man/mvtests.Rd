\name{mvtests}
\alias{mvtests}
\title{Multivariate SNP tests
}
\description{
  This function calculates multivariate score tests between a multivariate
  (or multinomial) phenotype and sets of SNPs
}
\usage{
mvtests(phenotype, sets, stratum, data = sys.parent(), snp.data, rules = NULL, complete = TRUE, uncertain = FALSE, score = FALSE)
}
\arguments{
  \item{phenotype}{
    Either a factor (for a multinomial phenotype) or a matrix (for a
    multivariate phenotype)
  }
  \item{sets}{
    A list of sets of SNPs to be tested against the phenotype
  }
  \item{stratum}{
    (Optional) a stratifying variable
  }
  \item{data}{
    A data frame in which \code{phenotype} and \code{stratum} reside. If
    absent these are assumed to be in the parent frame and correctly
    aligned with the rows of \code{snp.data}
  }
  \item{snp.data}{
    An object of class \code{SnpMatrix} containing the SNP data
  }
  \item{rules}{
    (Optional) A set of imputation rules. The function then carries out
    tess on imputed SNPs
  }
  \item{complete}{
    If \code{TRUE} each test will use only subjects who have complete
    data for the phenotype and all SNPs in the set to be tested. If
    \code{FALSE}, then complete data for the phenotype is required, but
    tests are based upon complete pairs of SNPs
  }
  \item{uncertain}{
    If \code{TRUE}, uncertain genotype calls will be used in the
    tests (scored by their posterior expectations).
    Otherwise such calls are treated as missing
  }
  \item{score}{
    If \code{TRUE}, the score vectors and their variance-covariance
    matrices are saved in the output object for further processing
  }
}
\details{
  Currently \code{complete=FALSE} is not implemented
}
\value{
  An object of class \code{\link[=GlmTests-class]{snp.tests.glm}}
  or \code{\link[=GlmTestsScore-class]{GlmTests.score}}
  depending on whether \code{score} is set to \code{FALSE} or \code{TRUE}
  in the call
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  This is an experimental version
}
\examples{
## No example yet
}
\keyword{htest}

