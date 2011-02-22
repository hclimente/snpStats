\name{GlmTests-class}
\docType{class}
\alias{GlmTests-class}
\alias{[,GlmTests,ANY,missing,missing-method}
\alias{chi.squared,GlmTests,missing-method}
\alias{deg.freedom,GlmTests-method}
\alias{effect.sign,GlmTests,logical-method}
\alias{names,GlmTests-method}
\alias{coerce,GlmTests,data.frame-method}
\alias{p.value,GlmTests,missing-method}
\alias{sample.size,GlmTests-method}
\alias{show,GlmTests-method}
\alias{summary,GlmTests-method}

\alias{GlmTestsScore-class}
\alias{[,GlmTestsScore,ANY,missing,missing-method}
\alias{pool2,GlmTestsScore,GlmTestsScore,logical-method}
\alias{switch.alleles,GlmTestsScore,character-method}
  
\title{Classes "GlmTests" and "GlmTestsScore"}
\description{Classes of objects created by \code{\link{snp.lhs.tests}}
  and  \code{\link{snp.rhs.tests}}. The class
    \code{"GlmTestsScore"} extends the class
    \code{"GlmTests"} and is invoked by setting the argument
    \code{score=TRUE} when calling testing functions
    in order to save the scores and their variances
    (and covariances)
}
\section{Objects from the Class}{Objects of class \code{"GlmTests"}
  have four slots:
  \describe{
    \item{snp.names}{When only single SNPs are tested,
      a  character vector of SNP names. Otherwise a list
      of such vectors (one for each test)}
    \item{var.names}{A character vector containing names of variables
      tested against SNPs}
    \item{chisq}{A numerical vector of chi-squared test values}
    \item{df}{An integer vector of degrees of freedom for the tests}
    \item{N}{A integer vector of the number of samples contributing to
      each test}
  }
  The \code{"GlmTestsScore"} class extends this, adding a slot
  \code{score} containing a list with elements which are themselves
  lists with two elements:
  \describe{
    \item{U}{The vector (or matrix) of efficient scores}
    \item{V}{The upper triangle of the variance-covariance matrix of
      \code{U}, stored as a vector}
  }
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x = "GlmTests", i = "ANY", j =
	"missing", drop = "missing")}: Subsetting operator}
    \item{coerce}{\code{signature(from = "GlmTests", to = "data.frame")}:
      Simplify object}
    \item{chi.squared}{\code{signature(x = "GlmTests", df =
	"missing")}: Extract chi-squared test values}
    \item{deg.freedom}{\code{signature(x = "GlmTests")}: Extract
      degrees of freedom for tests}
    \item{names}{\code{signature(x="GlmTests")}: Extract (or generate)
      a name for each test}
    \item{p.value}{\code{signature(x = "GlmTests", df =
	"missing")}: Extract \emph{p}-values}
    \item{sample.size}{\code{signature(object = "GlmTests")}:
      Extract sample sizes for tests}
    \item{show}{\code{signature(object = "GlmTests")}: Show method}
    \item{summary}{\code{signature(object = "GlmTests")}: Summary
      method}

    \item{[]}{\code{signature(x = "GlmTestsScore", i = "ANY", j =
	"missing", drop = "missing")}: Subsetting operator}
    \item{effect.sign}{\code{signature(x = "GlmTestsScore", simplify =
	"logical")}: Extract signs of associations. If \code{simpify} is
      \code{TRUE} then a simple vector is returned if all tests are on 1df}
    \item{pool2}{\code{signature(x = "GlmTestsScore", y =
	"GlmTestsScore", score = "logical")}: Combine results from two
      sets of tests}
    \item{switch.alleles}{\code{signature(x = "GlmTestsScore", snps =
	"character")}: Emulate, in the score vector and its (co)variances,
      the effect of switching of the alleles of specified SNPs}
  }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{Most of the methods for this class are shared with
the \code{\linkS4class{SingleSnpTests}} and
\code{\linkS4class{SingleSnpTestsScore}} classes
}
\seealso{\code{\link{snp.lhs.tests}},\code{\link{snp.rhs.tests}},
   \code{\linkS4class{SingleSnpTests}},
   \code{\linkS4class{SingleSnpTestsScore}}
}
\examples{
showClass("GlmTests")
}
\keyword{classes}
