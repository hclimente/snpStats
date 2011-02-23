\name{for.exercise}
\docType{data}
\alias{for.exercise}
\alias{snps.10}
\alias{snp.support}
\alias{subject.support}
\title{Data for exercise in use of the snpStats package}
\description{
These data have been created artificially from publicly available
datasets. The SNPs have been selected from those genotyped by the
International HapMap Project
(\url{http://www.hapmap.org})
to represent
the typical density found on a whole genome association chip, (the Affymetrix
500K platform,
\url{http://www.affymetrix.com/support/technical/sample\_data/500k\_hapmap\_genotype\_data.affx}
for a moderately sized chromosome (chromosome 10). A study of 500 cases and
500 controls has been simulated allowing for recombination using beta
software from Su and Marchini
(\url{http://www.stats.ox.ac.uk/~marchini/software/gwas/hapgen.html}).
Re-sampling of cases was weighted in such a
way as to simulate three ``causal'' locus on this chromosome, with
multiplicative effects of 1.3, 1.4 and 1.5 for each copy of the risk allele.
}
\usage{data(for.exercise)}
\format{
  There are three data objects in the dataset:
  \itemize{
    \item \code{snps.10}: 
      An object of class "\code{SnpMatrix}"
      containing a matrix of SNP genotype calls. Rows of the matrix
      correspond to subjects and columns correspond to SNPs.
    \item \code{snp.support}: 
      A conventional R
      data frame containing information about the
      SNPs typed (the chromosome position and the nucleotides
      corresponding to the two alleles of the SNP).
    \item \code{subject.support}:
      A conventional R dataframe containing information about the study
      subjects. There are two variables; \code{cc} gives case/control
      status (1=case), and \code{stratum} gives ethnicity.
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
data(for.exercise)
snps.10
summary(snps.10)
summary(snp.support)
summary(subject.support)
}
\keyword{datasets}
