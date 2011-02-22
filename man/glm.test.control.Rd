\name{glm.test.control}
\alias{glm.test.control}
\title{Set up control object for GLM tests}
\description{
 To carry out a score test for a GLM, we first fit a "base" model using
 the standard iteratively reweighted least squares (IRLS) algorithm and
 then carry out a score test for addition of further terms. This
 function sets various control parameters for this.
}
\usage{
glm.test.control(maxit, epsilon, R2Max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxit}{Maximum number of IRLS steps}
  \item{epsilon}{Convergence threshold for IRLS algorithm}
  \item{R2Max}{R-squared limit for aliasing of new terms}
}
\details{
Sometimes (although not always), an iterative scheme is necessary to fit
the "base"  generalized linear model (GLM) before carrying out a score
test for effect of adding new term(s). The \code{maxit} parameter sets
the maximum number of iterations to be carried out, while the
\code{epsilon} parameter sets the criterion for determining
convergence. After fitting the base model, the new terms are added, but
terms judged to be "aliased" are omitted. The method for determining
aliasing is as follows (denoting the "design" matrix for the additional
terms by \code{Z}):
  \enumerate{ 
    \item{Step 1}{Regress each column of \code{Z} on the base model matrix,
    using the final GLM weights from the base model fit, and replace
    \code{Z} with the residuals from these regressions.}
    \item{Step 2}{Consider each column of the new \code{Z} matrix in turn,
    regressing it on the \emph{previous} columns (again using the weights
    from the base model fit). If the proportion of the weighted sum of
    squares "explained" by this regression exceeds \code{R2Max}, the term
    is dropped and not included in the test,}
  }
The aim of this procedure to avoid wasting degrees of freedom on columns
so strongly aliased that there is little power to detect their effect.
}
\value{
  Returns the parameters as a list in the expected order
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.lhs.tests}}, \code{\link{snp.rhs.tests}}}
\keyword{utilities}
