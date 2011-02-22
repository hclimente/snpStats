\name{GlmEstimates-class}
\Rdversion{1.1}
\docType{class}
\alias{GlmEstimates-class}
\alias{[,GlmEstimates,ANY,missing,missing-method}
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
    \item{[}{\code{signature(x = "GlmEstimates", i = "ANY", j = "missing", drop = "missing")}: ... }
    \item{show}{\code{signature(object = "GlmEstimates")}: ... }
	 }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{
  \code{\link{snp.lhs.estimates}},
  \code{\link{snp.rhs.estimates}}
}
\examples{
showClass("GlmEstimates")
}
\keyword{classes}
