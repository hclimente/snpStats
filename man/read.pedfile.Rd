\name{read.pedfile}
\alias{read.pedfile}
\title{
Read a pedfile as \code{"SnpMatrix"} object
}
\description{
Reads diallelic data in linkage "pedfile" format, with one line of data
per sample (subject) containing six mandatory fields followed by pairs
of fields, one pair for each locus,  giving the two alleles observed.
}
\usage{
read.pedfile(file, n, snps, which, split = "\t| +", sep = ".", na.strings = "0", lex.order = FALSE)
}
\arguments{
  \item{file}{
    The input pedfile. This may be (but need not be) gzipped
}
  \item{n}{
    (Optional) The number of lines of data to be read. If not supplied
    the pedfile is read once and rewound to determine how many lines it
    contains 
}
  \item{snps}{
    (Optional) Either a character vector giving the names of the loci,
    or a single character variable giving the name of a locus
    information file from which these can be read. This file is assumed
    to be white-space delimited with one line per locus and no header
    line. If this argument is not supplied, locus names are generated as
    a numerical sequence, prefixed by \code{locus} and a separator character
}
  \item{which}{
    (Optional) If locus names are to be read from a file, this argument
    should specify which column contains the names. If not supplied, the
    first column giving unique locus names is used
}
  \item{split}{
    A "regexp" specifying how the input pedfile will be split into
    fields. The default value specifies either a TAB character or one or
    more spaces
}
  \item{sep}{
    The separator character used in constructing row and column names of
    the output \code{SnpMatrix} object
}
  \item{na.strings}{
    One or more strings to be set to \code{NA}. Any field taking one  of these
    values will be set to \code{NA}
  }
  \item{lex.order}{
    If TRUE, then alleles will be allocated to internal \code{1} and
    \code{2} values in lexographic order. Otherwise they are converted
    in the order in which they are encountered when reading the file
    (the default setting)
}
}
\details{
  Row names for the output \code{SnpMatrix} object and for the
  accompanying subject description dataframe are taken as the pedigree
  identifiers, when these provide the required unique identifiers. When
  these are duplicated, an attempt is made to use the pedigree-member
  identifiers instead but, when these too are duplicated,  
  row names are obtained by concatenating, with a separator character, the
  pedigree and pedigree-member identifiers.
}
\value{
A list, comprising
\item{genotypes}{The output genotype data as an object of class
  \code{"SnpMatrix"}. If either the pedigree or pedigree-member identifiers in
  the \code{ped} file are not duplicated, these are used for the row
  names of the output object. Otherwise these two fields are
  concatenated, separated by \code{sep}}
\item{fam}{A dataframe containing the first six fields in the
  pedfile. The row names will correspond with those of the \code{SnpMatrix}}
\item{map}{A dataframe giving the alleles at each locus. If locus names
  were obtained from a dataframe read from an existing file, then the
  allele information is simply appended to this frame. Otherwise a new
  dataframe is created. The row names will correspond with the column
  names of the \code{SnpMatrix }}
}
\note{
  This function is written entirely in R and may not be particularly
  fast. However, it imposes no restrictions on the allele codes
  recognized.

  Homozygous genotypes may be represented in the input file either (a)
  by coding both alleles to the same value, or (b) setting the second
  allele to "missing" (as specified by the \code{missing.allele}
  argument). No special provision is made to read \code{XSnpMatrix}
  objects; such data should first be read as a \code{SnpMatrix} and then
  coerced to an \code{XSnpMatrix} using \code{new} or \code{as}.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}}
\examples{
##
## No example supplied yet
##
}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
