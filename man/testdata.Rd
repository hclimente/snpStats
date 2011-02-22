\name{testdata}
\docType{data}
\alias{testdata}
\alias{Autosomes}
\alias{Xchromosome}
\alias{Asnps}
\alias{Xsnps}
\alias{subject.data}
\title{Test data for the snpAssoc package}
\description{
  This dataset comprises several data frames from a fictional (and
  unrealistically small) study. The
  dataset started off as real data from a screen of non-synonymous SNPs
  for association with type 1 diabetes, but the original identifiers have
  been removed and a random case/control status has been generated.
}
\usage{data(testdata)}
\format{
  There are five data objects in the dataset:
  \itemize{
    \item \code{Autosomes}: An object of class \code{"SnpMatrix"}
    containing genotype calls for 400 subjects at 9445 autosomal SNPs
    \item \code{Xchromosome}: An object of class \code{"XSnpMatrix"}
    containing genotype calls for 400 subjects at 155 SNPs on the X
    chromosome
    \item \code{Asnps}: A dataframe containing information about the
    autosomal SNPs. Here it contains only one variable,
    \code{chromosome}, indicating the chromosomes on which the SNPs are
    located
    \item \code{Xsnps}: A dataframe containing information about the X
    chromosome SNPs. Here it is empty and is only included for
    completeness
    \item \code{subject.data}: A dataframe containing information about
    the subjects from whom each row of SNP data was obtained. Here it
    contains:
    \itemize{
      \item \code{cc}: Case-control status
      \item \code{sex}: Sex
      \item \code{region}: Geographical region of residence
    }
  }
}
\source{
  The data were obtained from the diabetes and inflammation laboratory
  (see \url{http://www-gene.cimr.cam.ac.uk/todd})
}
\references{
  \url{http://www-gene.cimr.cam.ac.uk/clayton}
}
\examples{
data(testdata)
Autosomes
Xchromosome
summary(Asnps)
summary(Xsnps)
summary(subject.data)
summary(summary(Autosomes))
summary(summary(Xchromosome))
}
\keyword{datasets}
