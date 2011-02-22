\name{read.beagle}
\alias{read.beagle}
\title{
Read genotypes imputed by the BEAGLE program
}
\description{
The BEAGLE program generates, for each SNP and each subject, 
posterior probabilities for the three genotypes. This function reads
such data as a \code{SnpMatrix} object, storing the posterior
probabilities to as much accuracy allowed by a one-byte coding
}
\usage{
read.beagle(file, rownames=NULL, nsnp = NULL, header=TRUE)
}
\arguments{
  \item{file}{
    The input file name. This file my be gzipped.
  }
  \item{rownames}{
    The row names (sample identifiers) for the matrix
    }
  \item{nsnp}{
    The number of SNPs to be read in. This corresponds with the number
    of lines in the input file. If not supplied, the function does a
    preliminary pass to determine the number of lines
  }
  \item{header}{
    Set this \code{TRUE} if the file contains a header line (it won't
    for older versions of BEAGLE)
  }
}
\details{
  In later versions of BEAGLE, row names are listed on a header line.
  However, if the \code{rownames} argument is supplied, this will
  take precedence over the header line. If there is no header
  line and no row names are supplied, names
  are generated as \code{Sample1}, \code{Sample2} etc.
  
  No provision is made for data for the X chromosome. Such data must be
  first read as a \code{SnpMatrix} and subsequently coerced to an
\code{XSnpMatrix} object
}
\value{
an object of class \code{SnpMatrix}
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{SnpMatrix-class}}}
\examples{
##---- No example available yet
}
\keyword{IO}
\keyword{file}
