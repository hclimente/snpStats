\name{snp.lhs.tests}
\alias{snp.lhs.tests}
\title{Score tests with SNP genotypes as dependent variable}
\description{
  Under the assumption of Hardy-Weinberg equilibrium, a SNP genotype is
  a binomial variate with two trials for an autosomal SNP or with one or
  two trials (depending on sex) for a SNP on the X chromosome.
  With each SNP in an input
  \code{"SnpMatrix"} as dependent variable, this function first fits a
  "base" logistic regression model and then carries out a score test for
  the addition of further term(s). The Hardy-Weinberg
  assumption can be relaxed by use of a "robust" option.
}
\usage{
snp.lhs.tests(snp.data, base.formula, add.formula, subset, snp.subset,
                data = sys.parent(), robust = FALSE, uncertain = FALSE, 
                control=glm.test.control(maxit=20, epsilon=1.e-5,R2Max=0.999),
                score=FALSE)
}
\arguments{
  \item{snp.data}{The SNP data, as an object of class
    \code{"SnpMatrix"} or \code{"XSnpMatrix"} }
  \item{base.formula}{A \code{formula} object describing the base model,
    with dependent variable omitted }
  \item{add.formula}{A \code{formula} object describing the additional
    terms to be tested, also with dependent variable omitted}
  \item{subset}{An array describing the subset of observations to be
    considered} 
  \item{snp.subset}{An array describing the subset of SNPs to be
    considered. Default action is to test all SNPs.} 
  \item{data}{The data frame in which \code{base.formula},
    \code{add.formula} and \code{subset} are to be evaluated}
  \item{robust}{If \code{TRUE}, a test which does not assume
    Hardy-Weinberg equilibrium will be used }
  \item{uncertain}{If \code{TRUE}, uncertain genotypes are used and
    scored by their posterior expectations. Otherwise they are treated
    as missing. If set, this option forces \code{robust} variance estimates}
  \item{control}{An object giving parameters for the IRLS algorithm
    fitting of the base model and for the acceptable aliasing amongst
    new terms to be tested. See \code{\link{glm.test.control}}}
  \item{score}{Is extended score information to be returned?}
}
\details{
  The tests used are asymptotic chi-squared tests based on the vector of
  first and second derivatives of the log-likelihood with respect to the
  parameters of the additional model. The "robust" form is a generalized
  score test in the sense discussed by Boos(1992). 
  If a \code{data} argument is supplied, the \code{snp.data} and
  \code{data} objects are aligned by rowname. Otherwise all variables in
  the model formulae are assumed to be stored in the same order as the
  columns of the \code{snp.data} object. 
}
\value{
  An object of class \code{\link[=GlmTests-class]{snp.tests.glm}}
  or \code{\link[=GlmTestsScore-class]{GlmTests.score}}
  depending on whether \code{score} is set to \code{FALSE} or \code{TRUE}
  in the call.
}
\references{Boos, Dennis D. (1992) On generalized score tests. \emph{The
  American Statistician}, \strong{46}:327-333.}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  A factor (or
  several factors) may be included as arguments to the function
  \code{strata(...)} in the \code{base.formula}. This fits all
  interactions of the factors so included, but leads to faster
  computation  than fitting these in the normal way. Additionally, a
  \code{cluster(...)} call may be included in the base model
  formula. This identifies clusters of potentially correlated
  observations (e.g. for members of the same family); in this case, an
  appropriate robust estimate of the variance of the score test is used.
  No more than one
  \code{strata()} call may be used, and neither \code{strata(...)}  or
  \code{cluster(...)} calls may appear in the \code{add.formula}.
  A known bug is that the function fails when no \code{data} argument is
  supplied and the base model formula contains no variables
  (\code{~1}). A work-round is to create a data frame to hold the
  variables in the models and pass this as \code{data=}. 
}
\seealso{\code{\link{GlmTests-class}},
  \code{\link{GlmTestsScore-class}},
  \code{\link{glm.test.control}},\code{\link{snp.rhs.tests}}
    \code{\link{single.snp.tests}}, \code{\link{SnpMatrix-class}},
    \code{\link{XSnpMatrix-class}}}
\examples{
data(testdata)
snp.lhs.tests(Autosomes[,1:10], ~cc, ~region, data=subject.data)
snp.lhs.tests(Autosomes[,1:10], ~strata(region), ~cc,
   data=subject.data)
}
\keyword{htest}
