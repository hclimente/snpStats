\name{switch.alleles}
\alias{switch.alleles}
\title{Switch alleles in columns of a SnpMatrix or in test results}
\description{
  This is a generic function which can be applied to objects of class
  \code{"SnpMatrix"} or  \code{"XSnpMatrix"} (which hold SNP genotype
  data), or to objects of class 
  \code{"SingleSnpTestsScore"} or \code{"GlmTests"} (which hold
  association test results). In the former case, specified SNPs can be
  recoded as if the alleles were switched (so that \emph{AA} genotypes
  become \emph{BB} and 
  vice-versa while \emph{AB} remain unchanged). In the latter case,
  test results are modified \emph{as if} alleles had been switched. 
}
\usage{
switch.alleles(x, snps)
}
\arguments{
  \item{x}{The input object, of class  \code{"SnpMatrix"},
    \code{"XSnpMatrix"}, \code{"SingleSnpTestsScore"}, or
    \code{"GlmTests"}}
  \item{snps}{A vector of type integer, character or logical specifying
      the SNP to have its alleles switched}
}
\value{
  An object of the same class as the input object
}
\note{
  Switching alleles for SNPs has no effect on test results. These
  functions are required when carrying out meta-analysis, bringing
  together several sets of results. It is then important that alleles
  line up in the datasets to be combined. It is often more convenient
  (and faster) to apply this process to the test result objects rather
  than to the genotype data themselves.
  }
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{  \code{\link{SnpMatrix-class}},
  \code{\link{XSnpMatrix-class}},
  \code{\link{SingleSnpTests-class}},
  \code{\link{GlmTests-class}}
}
\examples{
data(testdata)
which <-  c("173774", "173811")
Asw <- switch.alleles(Autosomes, which)
col.summary(Autosomes[,which])
col.summary(Asw[,which])
}
\keyword{utilities}

