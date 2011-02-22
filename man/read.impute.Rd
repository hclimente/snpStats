\name{read.impute}
\alias{read.impute}
\title{
Read genotypes imputed by the IMPUTE2 program
}
\description{
The IMPUTE2 program generates, for each SNP and each subject, 
posterior probabilities for the three genotypes. This function reads
such data as a \code{SnpMatrix} object, storing the posterior
probabilities to as much accuracy allowed by a one-byte coding
}
\usage{
read.impute(file, rownames = NULL, nsnp = NULL, snpcol = 2)
}
\arguments{
  \item{file}{
    The input file name. This file my be gzipped.
  }
  \item{rownames}{
    The row names for the output object. Note that these correspond to
    groups of three columns in the input file. If not supplied, names
    are generated as \code{Sample1}, \code{Sample2} etc.
  }
  \item{nsnp}{
    The number of SNPs to be read in. This corresponds with the number
    of lines in the input file. If not supplied, the function does a
    preliminary pass to determine the number of lines
  }
  \item{snpcol}{
    Which column of the input will be used as the SNP name. Default is
    column 2 
  }
}
\details{
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
