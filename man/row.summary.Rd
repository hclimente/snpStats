\name{row.summary}
\alias{row.summary}
\alias{col.summary}
\title{Summarize rows or columns of a snp matrix}
\description{
  These function calculates summary statistics of each row or column of
  call rates and heterozygosity for each row of a an
  object of class \code{"SnpMatrix"} or \code{"XSnpMatrix"}
}
\usage{
row.summary(object)
col.summary(object, rules = NULL, uncertain = TRUE)
}
\arguments{
\item{object}{genotype data as a \code{\link{SnpMatrix-class}} or
  \code{\link{XSnpMatrix-class}} object}
\item{rules}{An object of class
  \code{"ImputationRules"}. If
  supplied, the rules coded in this object are used, together with
  the snp genotype data in
  \code{object}, to generate imputed SNPs. The column summary of these
  imputed data are then returned}
\item{uncertain}{If \code{TRUE} uncertain genotypes are used in
  calculation of allele and genotype frequencies (by scoring as
  posterior expectations). Otherwise, and for Hardy-Weinberg tests,
  they are ignored}
}
\value{
  \item{row.summary}{
    returns a data frame with rows corresponding to rows of the input object
    and with columns/elements:
    \itemize{
      \item Call.rate: Proportion of SNPs called
      \item Certain.calls: Proportion of called SNPs with certain calls
      \item Heterozygosity: Proportion of called SNPs which are
	  heterozygous
	}
	Uncertain calls are ignored for calculating the heterozygosity.
    	
	
  }
  \item{col.summary}{
    returns a data frame with rows corresponding to columns of the input object
    and with columns/elements:
    \itemize{
      \item Calls: The number of valid calls
      \item Call.rate: The proportion of genotypes called
      \item Certain.calls: Proportion of called SNPs with certain calls
      \item RAF: The "risk" allele (allele \code{B}) frequency
      \item MAF: The minor allele frequency
      \item P.AA: The frequency of homozygous genotype 1 (A/A)
      \item P.AB: The frequency of heterozygous genotype 2 (A/B)
      \item P.BB: The frequency of homozygous genotype 3 (B/B)
      \item z.HWE: A z-test for Hardy-Weinberg equilibrium
    }
    For objects of class \code{"XSnpMatrix"}, the following
    additional columns are returned:
    \itemize{
      \item P.AY: The frequency of allele A in males
      \item P.BY: The frequency of allele B in males
      \item Calls.female: The number of valid calls in females (only
	  these calls are used in the z-test for HWE)
      }
  }
}
\note{The current version of \code{row.summary}
  does not deal with the X chromosome
  differently, so that males are counted as homozygous.}
\author{David Clayton \email{dc208@cam.ac.uk}}
\examples{
data(testdata)
rs <- row.summary(Autosomes)
summary(rs)
cs <- col.summary(Autosomes)
summary(cs)
cs <- col.summary(Xchromosome)
summary(cs)
}
\keyword{utilities}
