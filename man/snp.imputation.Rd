\name{snp.imputation}
\alias{snp.imputation}
\title{Calculate imputation rules}
\description{
  Given two set of SNPs typed in the same subjects, this function
  calculates rules which can be used to impute one set
  from the other in a subsequent sample. The function can also calculate
  rules for imputing each SNP in a single dataset from other SNPs in the
  same dataset
}
\usage{
snp.imputation(X, Y, pos.X, pos.Y, phase=FALSE, try=50, stopping=c(0.95, 4, 0.05),
               use.hap=c(1.0, 0.0), em.cntrl=c(50,0.01,10,0.01), minA=5)
}
\arguments{
  \item{X}{An object of class \code{"SnpMatrix"} or
    \code{"XSnpMatrix"} containing observations
    of the SNPs to be used for imputation ("predictor SNPs")}
  \item{Y}{An object of same class as \code{X} containing observations
  of the SNPs to be imputed in a future sample ("target SNPs"). If this
  argument is missing, then target SNPs are also drawn from \code{X}}
  \item{pos.X}{The positions of the  predictor SNPs. Can be missing if
    there is no \code{Y} argument and the columns of \code{X} are in
    genome position order}
  \item{pos.Y}{The positions of the target SNPs. Only required when
    a \code{Y} argument is present}
  \item{phase}{See "Details" below}
  \item{try}{The number of potential predictor SNPs to be
    considered in the stepwise regression procedure around each target
    SNP . The nearest \code{try} predictor SNPs to each target SNP
    will be considered}
  \item{stopping}{Parameters of the stopping rule for the stepwise
    regression (see below)}
  \item{use.hap}{Parameters to control use of the haplotype imputation
    method (see below)}
  \item{em.cntrl}{Parameters to control test for convergence of EM
    algorithm for fitting phased haplotypes (see below)}
  \item{minA}{A minimum data quantity measure for estimating pairwise
    linkage disequilibrium (see below)}
}
\details{
  The routine first carries out a series of step-wise least-square
  regression analyses in
  which each Y SNP is regressed on the nearest \code{try} predictor (X)
  SNPs. If
  \code{phase} is \code{TRUE}, the regressions will be calculated at the
  chromosome (haplotype) level, variances being simply  \eqn{p(1-p)} and
  covariances estimated from the estimated two-locus haplotypes (this option is
  not yet implemented). Otherwise, the
  analysis is carried out at the genotype level based on  
  conventional variance and covariance estimates using the
  \code{"pairwise.complete.obs"} missing value treatment
  (see \code{\link{cov}}). New
  SNPs are added to the regression until either (a) the value of
  \eqn{R^2} exceeds the first parameter of \code{stopping}, (b) the
  number of "tag" SNPs has reached the maximum set in the second parameter of
  \code{stopping}, or (c) the change in \eqn{R^2} does not achieve the
  target set by the third parameter of \code{stopping}. If the third
  parameter of \code{stopping} is \code{NA}, this last test is replaced
  by a test for improvement in the Akaike information criterion
  (AIC).

  After choosing the set of "tag" SNPs in this way, a prediction
  rule is generated either by calculating phased haplotype frequencies,
  either (a) under a log-linear model for linkage disequilibrium with
  only first order association terms fitted, or (b) under the
  "saturated" model.  
  These methods do not differ if there is only
  one tag SNP but, otherwise,  choice between  methods is controlled
  by the  \code{use.hap} parameters. 
  If the prediction,  as measure by  \eqn{R^2} achieved with the
  log-linear smoothing model exceeds a
  threshold (the first parameter of \code{use.hap})
  then this method is used. Otherwise, if the gain in \eqn{R^2}
  achieved by using the second method exceeds the second parameter of
  \code{use.hap}, then the second method is used.
  Current experience is that, the log-linear method is rarely
  preferred with reasonable choices for \code{use.hap}, and imputation
  is much faster when the second  method only is considered. 
  The current default ensures that this second method is used,
  but the other possibility might be considered if imputing
  from very small samples; however this code is not extensively tested
  and should be regarded as experimental.

  The argument \code{em.cntrl} controls convergence
  testing for the EM algorithm for fitting haplotype frequencies and the
  IPF algorithm for fitting the log-linear model. The
  first parameter is the maximum number of EM iterations, and the second
  parameter is the threshold for the change in log likelihood
  below which the iteration is judged to have converged. The third and
  fourth parameters give the maximum number of IPF iterations and the
  convergence tolerance. There should be no need to change the default
  values. 

  All SNPs selected for imputation must have sufficient data for
  estimating pairwise linkage disequilibrium with each other and with
  the target SNP. The statistic chosen is based on the four-fold tables
  of two-locus haplotype frequencies. If the frequencies in such a table
  are labelled \eqn{a, b, c} and \eqn{d} then, if \eqn{ad>bc} then
  \eqn{t = min(a,d)} and, otherwise,  \eqn{t = min(b,c)}. The cell
  frequencies \eqn{t} must exceed \code{minA} for all pairwise
  comparisons.  
}
\value{
  An object of class
  \code{"ImputationRules"}.
}
\references{
  Chapman J.M., Cooper J.D., Todd J.A. and Clayton D.G. (2003)
  \emph{Human Heredity}, \bold{56}:18-31.
  
  Wallace, C. et al. (2010) \emph{Nature Genetics}, \bold{42}:68-71
}
\note{The \code{phase=TRUE} option is not yet implemented} 
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{ImputationRules-class}},
  \code{\link{imputation.maf}}, \code{\link{imputation.r2}}}
\examples{
# Remove 5 SNPs from a datset and derive imputation rules for them
data(for.exercise)
sel <- c(20, 1000, 2000, 3000, 5000)
to.impute <- snps.10[,sel]
impute.from <- snps.10[,-sel]
pos.to <- snp.support$position[sel]
pos.fr <- snp.support$position[-sel]
imp <- snp.imputation(impute.from, to.impute, pos.fr, pos.to)
}
\keyword{models}
\keyword{regression}
