\name{ld}
\alias{ld}
\title{
  Pairwise linkage disequilibrium measures
}
\description{
  This function calculates measures of linkage disequilibrium between
  pairs of SNPs. The two SNPs in each pair may both come from the same
  \code{SnpMatrix} object, or from two different \code{SnpMatrix}
  objects. Statistics which can be calculated are the
  log likelihood ratio, odds ratio, Yule's Q, covariance, D-prime,
  R-squared, and R.
}
\usage{
ld(x, y = NULL, depth = NULL, stats, symmetric = FALSE)
}
\arguments{
  \item{x}{
    An object of class \code{SnpMatrix} or \code{XSnpMatrix}
  }
  \item{y}{
    (Optional) Another object of the same class as \code{x}. If \code{y}
    is supplied, LD statistics are calculated between each column of
    \code{x} and each column of \code{y}. Otherwise, they are calculated
    between columns of \code{x}
  }
  \item{depth}{
    When \code{y} is not supplied, this parameter is mandatory and
    controls the maximum lag between columns of \code{x}
    considered. Thus, LD statistics are calculated between \code{x[,i]}
    and \code{x[,j]} only if \code{i} and \code{j} differ by no more
    than \code{depth}
  }
  \item{stats}{
    A character vector specifying the linkage disequilibrium measures to
    be calculated. This should contain one or more of the strings: \code{"LLR"},
    \code{"OR"}, \code{"Q"}, \code{"Covar"}, \code{"D.prime"},
    \code{"R.squared"}, ad \code{"R"}
  }
  \item{symmetric}{
    When no \code{y} argument is supplied this argument controls the format of
    the output band matrices. If \code{TRUE}, symmetric matrices are
    returned and, otherwise, an upper triangular matrices are returned
  }
}
\details{
For each pair of SNPs, phased haplotype frequencies are first estimated
by maximum likelihood using the method described by Clayton and Leung
(2007). The arrays of chosen LD statistics are then calculated and
returned, either as
band matrices (when \code{y} is not supplied), or as conventional
rectangular matrices (when \code{y} is supplied). Band matrices are
stored in compressed form as objects of class \code{dsCMatrix}
(symmetric) or \code{dgCMatrix} (upper triangular). These classes are
defined in the \code{"\link[Matrix:Matrix-class]{Matrix}"} package)
}
\value{
  If only one LD statistic is requested, the function returns either a
  matrix or a compressed band matrix. If more than one LD statistic is
  requested, a list of such objects is returned
}
\references{
  Clayton and Leung (2007) \emph{Human Heredity}, \bold{64}:45-51, 
  (this paper is included in package documentation)
}
\author{
David Clayton \email{dc208@cam.ac.uk}
}
\seealso{
\code{"\link[Matrix:Matrix-class]{Matrix-class}"}
}
\examples{
data(testdata)
ld1 <- ld(Autosomes[, 1:50], depth=10, stats=c("D.prime", "R.squared"))
ld2 <- ld(Autosomes[, 1:20], Autosomes[, 21:25], stats="R.squared") 
}
\keyword{misc}

