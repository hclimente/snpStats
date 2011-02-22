\name{qq.chisq}
\alias{qq.chisq}
\title{Quantile-quantile plot for chi-squared tests}
\description{
This function plots ranked observed chi-squared test statistics against
the corresponding expected order statistics. It also estimates an
inflation (or deflation) factor, lambda, by the ratio of the trimmed
means of observed and expected values. This is useful for
inspecting the results of whole-genome association studies for
overdispersion due to population substructure and other sources of bias
or confounding.
}
\usage{
qq.chisq(x, df=1, x.max, main="QQ plot", 
    sub=paste("Expected distribution: chi-squared (",df," df)", sep=""), 
    xlab="Expected", ylab="Observed",
    conc=c(0.025, 0.975), overdisp=FALSE, trim=0.5,  
    slope.one=FALSE, slope.lambda=FALSE, pvals=FALSE,
    thin=c(0.25,50), oor.pch=24, col.shade="gray", ...)
}
\arguments{
  \item{x}{A vector of observed chi-squared test values}
  \item{df}{The degreees of freedom for the tests}
  \item{x.max}{If present, truncate the observed value (Y) axis at
    \code{abs(x.max)}. If \code{x.max} is negative, the y-axis will
    extend to \code{abs(x.max)} even if the observed data do not}
  \item{main}{The main heading}
  \item{sub}{The subheading}
  \item{xlab}{x-axis label (default "Expected")}
  \item{ylab}{y-axis label (default "Observed")}
  \item{conc}{Lower and upper probability bounds for concentration band
    for the plot. Set this to \code{NA} to suppress this}
  \item{overdisp}{If \code{TRUE}, an overdispersion factor, lambda, will be 
    estimated and used in calculating concentration  band} 
  \item{trim}{Quantile point for trimmed mean calculations for estimation of 
    lambda. Default is to trim at the median}
  \item{slope.one}{Is a line of slope one to be superimpsed?}
  \item{slope.lambda}{Is a line of slope lambda 
    to be superimposed?}
  \item{pvals}{Are P-values to be indicated on an axis drawn
     on the right-hand side of the plot?}
  \item{thin}{A pair of numbers indicating how points will be thinned
    before plotting (see Details).
    If \code{NA}, no thinning will be carried out}
  \item{oor.pch}{Observed values greater than \code{x.max} are plotted
    at \code{x.max}. This argument sets the plotting symbol to be used
    for out-of-range observations}
  \item{col.shade}{The colour with which the concentration band will be
    filled}
  \item{...}{Further graphical parameter settings to be passed to
    \code{points()}}
}
\details{
  To reduce plotting time and the size of plot files, the smallest
  observed and expected points are thinned so that only a reduced number of
  (approximately equally spaced) points are plotted. The precise
  behaviour is controlled by the parameter
  \code{thin}, whose value should be a pair of numbers.
  The first number must lie
  between 0 and 1 and sets the proportion of the X axis over which 
  thinning is to be applied. The second number should be an integer and
  sets the maximum number of points to be plotted in this section.
  
  The "concentration band" for the plot is shown in grey. This region is 
  defined by upper and lower probability bounds for each order statistic. 
  The default is to use the 2.5% and 97.5% bounds.
  Note that this is not a simultaneous confidence region; the probability 
  that the plot will stray outside the band at some point exceeds 95%.

  When required, the dispersion factor is estimated by the ratio 
  of the observed trimmed mean to its expected value under the chi-squared 
  assumption.
}
\value{
  The function returns the number of tests, the number of values omitted
  from the plot (greater than \code{x.max}), and the estimated
  dispersion factor, lambda. 
}
\note{
  All tests must have the same number of degrees of freedom. If this is
  not the case, I suggest transforming to p-values and then plotting 
  -2log(p) as chi-squared on 2 df. 

  }
\references{Devlin, B. and Roeder, K. (1999) Genomic control for
  association studies. \emph{Biometrics}, \bold{55}:997-1004}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{single.snp.tests}}, \code{\link{snp.lhs.tests}},
  \code{\link{snp.rhs.tests}} }
\examples{
## See example the single.snp.tests() function
}
\keyword{hplot}

