\name{chi.squared}
\alias{chi.squared}
\alias{deg.freedom}
\alias{effect.sign}
\alias{p.value}
\alias{sample.size}
\alias{effective.sample.size}
\title{Extract test statistics and p-values}
\description{
Generic functions to extract values from the SNP
association test objects returned by various testing functions
}
\usage{
chi.squared(x, df)
deg.freedom(x)
effect.sign(x, simplify)
p.value(x, df)
sample.size(x)
effective.sample.size(x)
}
\arguments{
  \item{x}{An object of class \code{"SingleSnpTests"}, 
     \code{"SingleSnpTestsScore"}, or \code{"GlmTests"}}
   \item{df}{Either the numeric value 1 or 2 (not used when \code{x} is
     of class \code{"GlmTests"})}
   \item{simplify}{This switch is relevant when \code{x} is of class
     \code{"GlmTests"} and plays the same role as it does in
     \code{\link{sapply}}. If \code{simplify=TRUE}, where possible
     the output is returned as a simple numeric vector rather than as a
     list}
}
\details{
  These functions operate on objects created by 
  \code{\link{single.snp.tests}}, \code{\link{snp.lhs.tests}}, and
  \code{\link{snp.lhs.tests}}.

  The functions \code{chi.squared} and \code{p.value} return the
  chi-squared statistic and the corresponding \emph{p}-value.
  The argument \code{df} is only used for
  output from \code{single.snp.tests}, since this function calculates
  both 1 df and 2 df tests for each SNP. The functions
  \code{snp.lhs.tests} and \code{snp.rhs.tests} potentially calculate
  chi-squared tests on varying degrees of freedom, which can be
  extracted with \code{deg.freedom}. The function \code{effect.sign}
  indicates the direction of 
  associations. When applied to an output object from
  \code{snp.single.tests}, it returns \code{+1} if the association, as
  measured by the 1 df test, is positive and \code{-1} if the
  association is negative. Each test calculated by \code{GlmTests}
  are potentially tests of several parameters so that the effect sign
  can be a vector. Thus \code{effect.sign} returns a list of sign
  vectors unless, if \code{simplify=TRUE}, and it can be simplified as a
  single vector with one sign for each test.
  The function \code{sample.size} returns the number of observations
  actually used in the test, after exclusions due to missing data have
  been applied, and \code{effective.sample.size} returns the effective
  sample size which is less than the true sample size for tests on
  imperfectly imputed SNPs.
}
\value{
  A numeric vector containing the
  chi-squared test statistics or p-values. The output vector has a \code{names}
  attribute.
}
\note{The \code{df} and \code{simplify} arguments are not always
  required (or legal). See above}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{single.snp.tests}}, 
  \code{\link{snp.lhs.tests}}, \code{\link{snp.rhs.tests}},
  \code{\link{SingleSnpTests-class}},
  \code{\link{SingleSnpTestsScore-class}},
  \code{\link{GlmTests-class}}
}
\examples{
data(testdata)
tests <- single.snp.tests(cc, stratum=region, data=subject.data,
   snp.data=Autosomes, snp.subset=1:10)
chi.squared(tests, 1)
p.value(tests, 1)
}
\keyword{utilities}

