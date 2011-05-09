\name{snp.rhs.estimates}
\alias{snp.rhs.estimates}
\title{Fit GLMs with SNP genotypes as independent variable(s)}
\description{
This function fits a generalized linear model with phenotype as
dependent variable and with a series of SNPs (or small sets of SNPs) as
predictor variables.  Optionally, one or more potential confounders of
a phenotype-genotype association may be included in the model. 
In order to protect against misspecification
of the variance function,  "robust" estimates of the variance-covariance
matrix of estimates may be calculated in place of the usual
model-based estimates.
}
\usage{
snp.rhs.estimates(formula, family = "binomial", link, weights, subset, data
= parent.frame(), snp.data, sets=NULL, robust = FALSE,
   control=glm.test.control(maxit=20, epsilon=1.e-4, R2Max=0.98))
}
\arguments{
  \item{formula}{The model formula, with phenotype as dependent variable
  and any potential confounders as independent variables. Note that
  parameter estimates are not returned for these model terms}
  \item{family}{A string defining the generalized linear model
    family. This currently should (partially)  match one of
    \code{"binomial"}, \code{"Poisson"}, \code{"Gaussian"} or
      \code{"gamma"} (case-insensitive)}
  \item{link}{A string defining the link function for the GLM. This
  currently should (partially) match one of \code{"logit"},
  \code{"log"}, \code{"identity"} or \code{"inverse"}. The
  default action is to use the "canonical" link for the family selected}
  \item{data}{The dataframe in which the model formula is to be interpreted}
  \item{snp.data}{An object of class \code{"SnpMatrix"} or
    \code{"XSnpMatrix"} containing the SNP data} 
  \item{sets}{Either a vector of SNP names (or  numbers) for the SNPs
    to be added to the model formula, or a logical vector of length
    equal to the number of columns in \code{snp.data}
    or a list of short vectors
    defining sets of SNPs to be included (see \code{Details})}
  \item{weights}{"Prior" weights in the generalized linear model}
  \item{subset}{Array defining the subset of rows of  \code{data} to use}
  \item{robust}{If \code{TRUE}, robust tests will be carried out}
  \item{control}{An object giving parameters for the IRLS algorithm
    fitting of the base model and for the acceptable aliasing amongst
    new terms to be tested. See \code{\link{glm.test.control}}}
}
\details{
  Homozygous SNP genotypes are coded 0
  or 2 and heterozygous genotypes are coded 1. For SNPs on the X
  chromosome, males are coded as homozygous females. For X SNPs, it will
  often be
  appropriate to include sex of subject in the base model (this is
  not done automatically).
  The "robust" option causes Huber-White estimates of the
  variance-covariance matrix of the parameter estimates to be
  returned. These protect against mis-specification of the variance
  function in the GLM, for example if binary or count data are
  overdispersed,
    
  If a \code{data} argument is supplied, the \code{snp.data} and
  \code{data} objects are aligned by rowname. Otherwise all variables in
  the model formulae are assumed to be stored in the same order as the
  columns of the \code{snp.data} object.

  Usually SNPs to be fitted in models will be referenced by name. However,
  they can
  also be referenced by number, indicating the
  appropriate column in the input \code{snp.data}.  They can also
  be referenced by a logical selection vector of length equal to the
  number of columns in \code{snp.data}.
}
\value{
  An object of class \code{\link[=GlmEstimates-class]{snp.estimates.glm}}
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  A factor (or
  several factors) may be included as arguments to the function
  \code{strata(...)} in the \code{formula}. This fits all
  interactions of the factors so included, but leads to faster
  computation  than fitting these in the normal way. Additionally, a
  \code{cluster(...)} call may be included in the base model
  formula. This identifies clusters of potentially correlated
  observations (e.g. for members of the same family); in this case, an
  appropriate robust estimate of the variance of the score test is used.
}
\seealso{\code{\link{GlmEstimates-class}},
    \code{\link{snp.lhs.estimates}},
    \code{\link{snp.rhs.tests}}, 
    \code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}}
\examples{
data(testdata)
test <- snp.rhs.estimates(cc~strata(region), family="binomial",
   data=subject.data, snp.data= Autosomes, sets=1:10)
print(test)
test2 <- snp.rhs.estimates(cc~region+sex, family="binomial",
   data=subject.data, snp.data= Autosomes, sets=1:10)
print(test2)
test.robust <- snp.rhs.estimates(cc~strata(region), family="binomial",
   data=subject.data, snp.data= Autosomes, sets=1:10, robust=TRUE)
print(test.robust)
}
\keyword{htest}
