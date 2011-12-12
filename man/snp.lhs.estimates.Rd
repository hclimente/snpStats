\name{snp.lhs.estimates}
\alias{snp.lhs.estimates}
\title{Logistic regression with SNP genotypes as dependent variable}
\description{
  Under the assumption of Hardy-Weinberg equilibrium, a SNP genotype is
  a binomial variate with two trials for an autosomal SNP or with one or
  two trials (depending on sex) for a SNP on the X chromosome.
  With each SNP in an input
  \code{"SnpMatrix"} as dependent variable, this function fits a
  logistic regression model. The Hardy-Weinberg
  assumption can be relaxed by use of a "robust" option.
}
\usage{
snp.lhs.estimates(snp.data, base.formula, add.formula, subset, snp.subset,
                data = sys.parent(), robust = FALSE, uncertain = FALSE,
                control=glm.test.control(maxit=20, epsilon=1.e-5,R2Max=0.999))
}
\arguments{
  \item{snp.data}{The SNP data, as an object of class
    \code{"SnpMatrix"} or \code{"XSnpMatrix"} }
  \item{base.formula}{A \code{formula} object describing a base model
    containing those terms which are to be fitted but for which
    parameter estimates are not required
    (the dependent variable is omitted from the model formula) }
  \item{add.formula}{A \code{formula} object describing the additional
    terms in the model for which parameter estimates are required (again,
    the dependent variable is omitted)}
  \item{subset}{An array describing the subset of observations to be
    considered} 
  \item{snp.subset}{An array describing the subset of SNPs to be
    considered. Default action is to test all SNPs.} 
  \item{data}{The data frame in which \code{base.formula},
    \code{add.formula} and \code{subset} are to be evaluated}
  \item{robust}{If \code{TRUE}, Hardy-Weinberg equilibrium will is not
    assumed in calculating the variance-covariance matrix of parameter
    estimates}
  \item{uncertain}{If \code{TRUE}, uncertain genotypes are used and
    scored by their posterior expectations. Otherwise they are treated
    as missing. If set, this option forces \code{robust} variance estimates}
  \item{control}{An object giving parameters for the IRLS algorithm
    fitting of the base model and for the acceptable aliasing amongst
    new terms to be tested. See \code{\link{glm.test.control}}}
}
\details{
  The model fitted is the union of the \code{base.formula} and
  \code{add.formula} models, although parameter estimates (and their
  variance-covariance matrix) are only
  generated for the parameters of the latter.
  The "robust" option causes a Huber-White "sandwich" estimate of the
  variance-covariance matrix to be used in place of the usual inverse
  second derivative matrix of the log-likelihood (which assumes
  Hardy-Weinberg equilibrium).
  If a \code{data} argument is supplied, the \code{snp.data} and
  \code{data} objects are aligned by rowname. Otherwise all variables in
  the model formulae are assumed to be stored in the same order as the
  columns of the \code{snp.data} object. 
}
\value{
  An object of class \code{\link[=GlmEstimates-class]{GlmEstimates}}
}
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
  appropriate robust estimate of the variance-covariance matrix of
  parameter estimates is calculated. No more than one
  \code{strata()} call may be used, and neither \code{strata(...)}  or
  \code{cluster(...)} calls may appear in the \code{add.formula}.

  If uncertain genotypes (e.g. as a result of imputation) are used, the
  interpretation of the regression coefficients is questionable.

  A known bug is that the function fails when no \code{data} argument is
  supplied and the base model formula contains no variables
  (\code{~1}). A work-round is to create a data frame to hold the
  variables in the models and pass this as \code{data=}. 
}
\seealso{\code{\link{GlmEstimates-class}}, \code{\link{snp.lhs.tests}}}
\examples{
data(testdata)
test1 <-
snp.lhs.estimates(Autosomes[,1:10], ~cc, ~region, data=subject.data)
test2 <-
snp.lhs.estimates(Autosomes[,1:10], ~strata(region), ~cc,
   data=subject.data)
test3 <-
snp.lhs.estimates(Autosomes[,1:10], ~cc, ~region, data=subject.data, robust=TRUE)
test4 <-
snp.lhs.estimates(Autosomes[,1:10], ~strata(region), ~cc,
   data=subject.data, robust=TRUE)
test5 <- snp.lhs.estimates(Autosomes[,1:10], ~region+sex, ~cc, data=subject.data, robust=TRUE)
print(test1)
print(test2)
print(test3)
print(test4)
print(test5)
}
\keyword{htest}
