\name{sm.compare}
\alias{sm.compare}
\title{
Compare two SnpMatrix objects
}
\description{
For quality control purposes, it is sometimes necessary to compare
genotype data derived from different sources. This function facilitates
this. 
}
\usage{
sm.compare(obj1, obj2, row.wise = TRUE, col.wise = TRUE)
}
\arguments{
  \item{obj1}{
    The first of the two \code{SnpMatrix} objects to be compared
  }
  \item{obj2}{
    The second \code{SnpMatrix} object
  }
  \item{row.wise}{
    Calculate comparison statistics aggregated in a row-wise manner
  }
  \item{col.wise}{
    Calculate column-wise comparison statistics
  }
}
\details{
  Initially row and column names of the two objects are compared to
  identify subsets of subjects and SNPs which they have in
  common. Then, every instance of a SNP genotype in the two objects are
  compared and agreements and disagreements counted by row and/or by
  column. 
}
\value{
  If only one of the row-wise and column-wise summaries are to be
  calculated, the return value is a matrix with rows defined by subjects
  or SNPs and columns giving counts of:
  \item{Agree}{Agreements (all)}
  \item{Disagree}{Disgreements (all)}
  \item{NA.agree}{Genotype coded \code{NA} in both objects}
  \item{NA.disagree}{Genotype coded \code{NA} in only one object}
  \item{Hom.agree}{Objects agree and genotype is homozygous}
  \item{Hom.switch}{Genotype coded as homozygous in both objects,
    but alleles switched}
  \item{Het.agree}{Genotype coded as heterozygous in both objects}
  \item{Het.Hom}{Genotype coded as heterozygous in one object and
    homozygous in the other}

  
  If both row-wise and column-wise summaries are computed (the default
  behaviour) , the function
  returns a list containing two matrices of the form described
  above. These are named \code{row.wise} and \code{col.wise}
}
\author{
David Clayton \email{david.clayton@cimr.cam.ac.uk}
}
\note{
  No special provision is yet made for objects of class
  \code{XSnpMatrix}, in which haploid calls are coded as homozygous.
}
\seealso{
  \code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}
}
\examples{
##
## No example yet available
##
}
\keyword{utilities}
