\name{misinherits}
\alias{misinherits}
\title{Find non-Mendelian inheritances in family data}
\description{
For SNP data in families, this function locates all subjects whose
parents are in the dataset and tests each SNP for non-Mendelian
inheritances in these trios.
}
\usage{
misinherits(ped, id, father, mother, data = sys.parent(), snp.data)
}
\arguments{
  \item{ped}{Pedigree identifiers}
  \item{id}{Subject identifiers}
  \item{father}{Identifiers for subjects' fathers}
  \item{mother}{Identifiers for subjects' mothers}
  \item{data}{A data frame in which to evaluate the previous four arguments}
  \item{snp.data}{An object of class \code{"SnpMatrix"} containing the SNP
    genotypes to be tested}
}
\details{
    The first four arguments are usually derived from  a "pedfile". If a 
  data frame is supplied for the \code{data} argument, the first four
  arguments will be evaluated in this frame. Otherwise they will be evaluated 
  in the calling environment. If the arguments are missing, they will be 
  assumed to be in their usual positions in the pedfile data frame 
  i.e. in columns one to four.  If the pedfile data are obtained from
  a dataframe, the row names of the \code{data} and \code{snp.data}
  files will be used to align the pedfile and SNP data. Otherwise, these
  vectors will be assumed to be in the same order as the rows of
  \code{snp.data}. 
}
\value{
  A logical matrix. Rows are subjects with any non-Mendelian
  inheritances and columns are SNPs with any non-Mendelian
  inheritances. The body of the matrix details whether each subject has
  non-Mendelian inheritance at each SNP. If a subject has no recorded
  genotype for a specific SNP, the corresponding element of the output
  matrix is set to NA.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{tdt.snp}}}
\examples{
data(families)
misinherits(data=pedData, snp.data=genotypes)
}
\keyword{manip}

