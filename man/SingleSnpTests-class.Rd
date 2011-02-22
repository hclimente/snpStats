\name{SingleSnpTests-class}
\docType{class}
\alias{SingleSnpTests-class}
\alias{[,SingleSnpTests,ANY,missing,missing-method}
\alias{chi.squared,SingleSnpTests,numeric-method}
\alias{names,SingleSnpTests-method}
\alias{p.value,SingleSnpTests,numeric-method}
\alias{sample.size,SingleSnpTests-method}
\alias{effective.sample.size,SingleSnpTests-method}
\alias{show,SingleSnpTests-method}
\alias{summary,SingleSnpTests-method}
\alias{coerce,SingleSnpTests,data.frame-method}
\alias{SingleSnpTestsScore-class}
\alias{[,SingleSnpTestsScore,ANY,missing,missing-method}
\alias{effect.sign,SingleSnpTestsScore,missing-method}
\alias{pool2,SingleSnpTestsScore,SingleSnpTestsScore,logical-method}
\alias{switch.alleles,SingleSnpTestsScore,ANY-method}


\title{Classes "SingleSnpTests" and "SingleSnpTestsScore"}
\description{These are classes to hold the objects created by
\code{\link{single.snp.tests}} and provide methods for extracting key
elements. The class \code{"SingleSnpTestsScore"} extends class
\code{"SingleSnpTests"}
to  include the score and score variance statistics in order to
provide methods for pooling results from several studies or parts of a study
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("SingleSnpTests", ...)} and
  \code{new("SingleSnpTestsScore", ...)} but, more usually, will be created
  by calls to \code{\link{single.snp.tests}}
}
\section{Slots}{
  \describe{
    \item{\code{snp.names}:}{The names of the SNPs tested, as they
      appear as column names in the original \code{SnpMatrix}}
    \item{\code{chisq}:}{A two-column
      matrix holding the 1 and 2 df association tests}
    \item{\code{N}:}{The numbers of
      observations included in each test}
    \item{\code{N.r2}:}{For tests on imputed SNPs, the product of \code{N} and
	the imputation \eqn{r^2}. Otherwise a zero-length object}
    \item{\code{U}:}{(class \code{"SingleSnpTestsScore"}) Score statistics}
    \item{\code{V}:}{(class \code{"SingleSnpTestsScore"}) Score variances}
  }
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x =  "SingleSnpTests", i = "ANY")}:
      Subsetting operator}
    \item{[]}{\code{signature(x =  "SingleSnpTestsScore", i = "ANY")}:
      Subsetting operator}
    \item{chi.squared}{\code{signature(x = "SingleSnpTests", df =
	"numeric")}:
      Extract 1- and 2-df chi-squared test values}
    \item{effect.sign}{\code{signature(x = "SingleSnpTestsScore", simplify =
	"missing")}: Extract signs of associations tested by the 1df tests}
     \item{names}{\code{signature(x="SingleSnpTests")}: Extract names of
      test values (\code{snp.names} slot)}
    \item{p.value}{\code{signature(x = "SingleSnpTests", df =
	"numeric")}:
      Evaluate 1- and 2-df test p-values}
    \item{show}{\code{signature(object = "SingleSnpTests")}:
      List all tests and p-values}
    \item{coerce}{\code{signature(from = "SingleSnpTests", to = "data.frame")}:
      Conversion to data frame class}
    \item{sample.size}{\code{signature(object = "SingleSnpTests")}:
      Extract sample sizes for tests}
    \item{effective.sample.size}{\code{signature(object = "SingleSnpTests")}:
      Extract effective sample sizes for tests. For imputed tests, these
      are the real sample sizes multiplied by the corresponding
      R-squared values for imputation}
    \item{summary}{\code{signature(object = "SingleSnpTests")}:
      Summarize all tests and p-values}
    \item{pool2}{\code{signature(x = "SingleSnpTestsScore",
	y = "SingleSnpTestsScore", score = "logical")}: Combine two
      sets of test results. Used recursively by {\code{\link{pool}}}}
    \item{switch.alleles}{\code{signature(x = "SingleSnpTestsScore", snps =
	"ANY")}: Emulate, in the score vector and its (co)variances,
      the effect of switching of the alleles for the specified tests}
  }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{
\code{\link{single.snp.tests}}, \code{\link{pool}}
}
\examples{
showClass("SingleSnpTests")
showClass("SingleSnpTestsScore")
}
\keyword{classes}
