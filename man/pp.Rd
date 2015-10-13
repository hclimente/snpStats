\name{pp}
\alias{pp}
\title{
Unpack posterior probabilities from one-byte codes
}
\description{
In \code{snpStats}, the three "posterior" probabilities corresponding to
the possible values of an uncertain genotype are packed into a single
byte code (with, of course, some loss in accuracy). This function, which
is provided as an aid to writing new functions, unpacks the posterior
probabilities from the single byte codes.  
}
\usage{
pp(x, transpose = FALSE)
}
\arguments{
  \item{x}{
    A vector, length N, which can be coerced into type \code{raw}
  }
  \item{transpose}{
    If \code{FALSE}, the result is an Nx3 matrix of posterior
    probabilities. If \code{TRUE}, a 3xN matrix is returned.
  }
}
\value{
  A numeric matrix
}
\author{
David Clayton \email{dc208@cam.ac.uk}
}
\examples{
##
## Read imputed data from a file produced by MACH 
##
path <- system.file("extdata/mach1.out.mlprob.gz", package="snpStats")
mach <- read.mach(path)
pp(mach[1:50, 10])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
