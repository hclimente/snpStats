\name{filter.rules}
\alias{filter.rules}
\title{Filter a set of imputation rules}
\description{
Determine which imputation rules are broken by removal of some SNPs from
a study. This function is needed because, when if it emerges that
genotyping of some SNPs is not reliable, necessitating their removal from
study, we would also wish to remove any SNPs imputed on the basis of
these unreliable SNPs.
}
\usage{
filter.rules(rules, snps.excluded, exclusions = TRUE)
}
\arguments{
  \item{rules}{An object of class \code{"ImputationRules"} containing
    a set of imputation rules}
  \item{snps.excluded}{The names of the SNPs whose
    removal is to be investigated}
  \item{exclusions}{If \code{TRUE}, the names of the imputed SNPs which
    would be lost by removal of the SNPs listed in
    \code{snps.excluded}. If \code{FALSE}, the names of the imputed SNPs
  which would \emph{not} be lost are returned}
}
\value{
  A character vector containing the names of imputed SNPs to be removed
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{ImputationRules-class}},
    \code{\link{snp.imputation}}
}
\examples{
# No example yet
}
\keyword{models}
\keyword{regression}

