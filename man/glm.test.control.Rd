\name{glm.test.control}
\alias{glm.test.control}
\title{Set up control object for GLM computations}
\description{
  Several commands depend on fitting a generalized linear model (GLM),
  using the standard 
  iteratively reweighted least squares (IRLS) algorithm.  This
  function sets various control parameters for this.
}
\usage{
glm.test.control(maxit = 20, epsilon = 1.e-5, R2Max = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxit}{Maximum number of IRLS steps}
  \item{epsilon}{Convergence threshold for IRLS algorithm}
  \item{R2Max}{R-squared limit for aliasing of new terms}
}
\details{
Sometimes (although not always), an iterative scheme is necessary to fit
a generalized linear model (GLM). The \code{maxit} parameter sets
the maximum number of iterations to be carried out, while the
\code{epsilon} parameter sets the criterion for determining
convergence. Variables which are judged to be "aliased" are dropped.
A variable is judged to be aliased if RSS/TSS is less than (1-R2Max),
where
\itemize{
  \item{RSS }{is the residual (weighted) sum of squares from the
    regression of that variable on the variables which precede it in the
    model formula (and any stratification defined in a strata()
    call in th emodel formula), and}
  \item{TSS }{is the total (weighted) sum of squared deviations of this
    variable from its mean (or, when a strata() call is present,
    from its stratum-specific  means).}
}
The weights used in this calculation
are the "working" weights of the IRLS algorithm.
}
\value{
  Returns the parameters as a list in the expected order
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{snp.lhs.tests}}, \code{\link{snp.rhs.tests}}}
\keyword{utilities}
