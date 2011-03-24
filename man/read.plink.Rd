\name{read.plink}
\alias{read.plink}
\title{Read a PLINK binary data file as a SnpMatrix}
\description{
The package PLINK saves genome-wide association data in groups of three
files, with the extensions \code{.bed}, \code{.bim}, and
\code{.fam}. This function reads these files and creates an object of
class \code{"SnpMatrix"}
}
\usage{
read.plink(bed, bim, fam, na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{The name of the
    file containing the packed binary SNP genotype data. It should have
    the extension \code{.bed}; if it doesn't, then this extension will
    be appended}
  \item{bim}{The file containing the SNP descriptions}
  \item{fam}{The file containing subject (and, possibly, family)
    identifiers. This is basically a tab-delimited "pedfile"}
  \item{na.strings}{Strings in .bam and .fam files to be recoded as NA}
  \item{sep}{A separator character for constructing unique subject
    identifiers}
  \item{select.subjects}{A numeric vector indicating a
    subset of subjects to be selected from the input file (see details)} 
  \item{select.snps}{Either a numeric or a character vector indicating a
   subset of SNPs to be selected from the input file (see details)}
}
\details{
  If the \code{bed} argument does not contain a filename with the file
  extension \code{.bed}, then this extension is appended to the
  argument. The remaining two arguments are optional; their default
  values are obtained by replacing the \code{.bed} filename extension by
  \code{.bim} and \code{.fam} respectively. See the PLINK documentation
  for the detailed specification of these files.

  The \code{select.subjects} or \code{select.snps} argument can be used
  to read a subset of the data. Use of \code{select.snps} requires that
  the \code{.bed} file is in SNP-major order (the default in
  PLINK). Likewise, use of \code{select.snps} requires that
  the \code{.bed} file is in individual-major order. Subjects are
  selected by their numeric order in the PLINK files, while SNPs are
  selected either by order or by name.  Note that
  the order of selected SNPs/subjects in the output objects
  will be the same as
  their order in the PLINK files. 

  Row names for the output \code{SnpMatrix} object and for the
  accompanying subject description dataframe are taken as the pedigree
  identifiers, when these provide the required unique identifiers. When
  these are duplicated, an attempt is made to use the pedigree-member
  identifiers instead but, when these too are duplicated,  
  row names are obtained by concatenating, with a separator character, the
  pedigree and pedigree-member identifiers.
}
\value{
  A list with three elements:
  \item{genotypes}{The output genotype data as an object of class}
  \code{"SnpMatrix"}. 
  \item{fam}{A dataframe corresponding to the \code{.fam} file,
  containing the first six fields in a standard pedfile.
  The row names will correspond with those of the \code{SnpMatrix}}
  \item{map}{A dataframe correponding to the \code{.bim} file. the row
    names correpond with the column names of the  \code{SnpMatrix}}
}
\note{
  No special provision is made to read \code{XSnpMatrix}
  objects; such data should first be read as a \code{SnpMatrix} and then
  coerced to an \code{XSnpMatrix} using \code{new} or \code{as}.
}
\references{PLINK: Whole genome association analysis toolset.
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/}
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{write.plink}},
  \code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
