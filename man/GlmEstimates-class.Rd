\name{GlmEstimates-class}
\Rdversion{1.1}
\docType{class}
\alias{GlmEstimates-class}
\alias{[,GlmEstimates,ANY,missing,missing-method}
\alias{coerce,GlmEstimates,GlmTests-method}
\alias{show,GlmEstimates-method}

\title{Class "GlmEstimates"}
\description{A simple class to hold output from
  \code{\link{snp.lhs.estimates}} and
  \code{\link{snp.rhs.estimates}}. Its main purpose is to provide a
  \code{show} method}
\section{Objects from the Class}{
Objects from this class are simple lists. Each element of the list is a
list giving the results of a generalized linear model fit, with elements:
\describe{
  \item{Y.var}{Name of the Y variable}
  \item{beta}{The vector or parameter estimates (with their names)}
  \item{Var.beta}{The upper triangle of the variance-covariance matrix
    of estimates, stored as a simple vector}
  \item{N}{The number of "units" used in the model fit}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "GlmEstimates", i = "ANY", j = "missing", drop = "missing")}: Subset }
    \item{coerce}{\code{signature(from = "GlmEstimates", to =
	"GlmTests")}: Calculate Wald tests }
    \item{show}{\code{signature(object = "GlmEstimates")}: Display}
	 }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{Wald tests calculated by coercing an object of class
  \code{GlmEstimates} to class \code{GlmTests} are \emph{asymptotically}
  equivalent to the "score" tests calculated by
  \code{\link{snp.lhs.tests}} and \code{\link{snp.rhs.tests}}, and both types
  of test are \emph{asymptotically} correct. But the asymptotic properties of
  Wald tests are not as good as those of score tests and, in
  circumstances such as low sample size or low minor allele frequency,
  Wald and score tests may differ substantially. In general score tests
  are to be preferred, but Wald tests are provided in \code{snpStats}
  because they are widely used.
}
\seealso{
  \code{\link{snp.lhs.estimates}},
  \code{\link{snp.rhs.estimates}}
}
\examples{
showClass("GlmEstimates")
}
\keyword{classes}
