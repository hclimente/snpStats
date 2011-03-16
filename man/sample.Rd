\name{sample.ped.gz}
\docType{data}
\alias{sample.info}
\alias{sample.bed}
\alias{sample.bim}
\alias{sample.fam}
\alias{sample-long.gz}
\alias{mach1.out.mlprob.gz}
\title{Sample datasets to illustrate data input}
\description{
  The first five files concern data on 20 diallelic loci on 120
  subjects. These data are 
  distributed with the Haploview package (Barrett et al., 2003). The
  sixth file contains a additional dataset of 18 SNPs in 100 subjects,
  coded in long format. These six files are used in the data input
  vignette. The final file is a sample imputed genotype
  dataset distributed with the MACH imputation package, and used in the
  imputation vignette.

  These files are stored in the \code{inst/extdata} relative to the
  package base. Full file names can be obtained using the
  \code{system.file} function.
}
\format{
  There are five files:
  \itemize{
    \item \code{sample.ped.gz}: A gzipped pedfile
    \item \code{sample.info}: An accompanying locus information file
    \item \code{sample.bed}: The corresponding PLINK \code{.bed} file
    \item \code{sample.bim}: The PLINK \code{.bim} file
    \item \code{sample.fam}: The PLINK \code{.fam} file
    \item \code{sample-long.gz}: A sample of long-formatted data
    \item \code{mach1.out.mlprob.gz}: An \code{mlprob} output file from
    the MACH genotype imputation program. This file contains, for each
    imputed genotype call, posterior probabilities for the three possible
    genotypes
  }
}
\source{
  \url{http://www.broadinstitute.org/scientific-community/science/programs/medical-and-population-genetics/haploview/downloads}
  \url{http://www.sph.umich.edu/csg/abecasis/MACH/download}
}
\references{
Barrett JC, Fry B, Maller J, Daly MJ.(2005)  Haploview: analysis
and visualization of LD and haplotype maps. \emph{Bioinformatics}, 2005 Jan 15,
[PubMed ID: 15297300]
}
\keyword{datasets}
