\name{read.long}
\alias{read.long}
\title{
Read SNP genotype data in long format
}
\description{
This function reads SNP genotype data from a file in which each line
refers to a single genotype call. Replaces the earlier function
\code{read.snps.long}.
}
\usage{
read.long(file, samples, snps,
            fields = c(snp = 1, sample = 2, genotype = 3, confidence = 4,
                       allele.A = NA, allele.B = NA),
            split = "\t| +", gcodes, no.call = "", threshold = NULL,
            lex.order = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
    Name(s) of file(s) to be read (can be gzipped)
}
  \item{samples}{
    Either a vector of sample identifiers, or the number of samples to
    be read. If a single file is to be read and this argument is
    omitted, the file will be scanned initially and all
    samples will be included
}
  \item{snps}{
    Either a vector of SNP identifiers, or the number of SNPs to
    be read. If a single file is to be read and this argument is
    omitted, the file will be scanned initially and all
    SNPs will be included
  }
  \item{fields}{
    A named vector giving the locations of the required fields. See
    Details below
  }
  \item{split}{
    A regular expression specifying how the input line will be split
    into fields. The default value specifies separation of fields by a
    TAB character, or by one or more blanks
}
  \item{gcodes}{
    When the genotype is read as a single field, this argument specifies
    how it is handled. See Details below.
}
  \item{no.call}{
    The string which indicates "no call" for either a genotype or (when
    the genotype is read as two allele fields) an allele
}
  \item{threshold}{
    A vector of length 2 giving the lower and higher acceptable limits
    for the confidence score
}
  \item{lex.order}{
    If \code{TRUE}, the alleles at each locus will be in lexographical
    order. Otherwise, ordering of alleles is arbitrary, depending on
    the order in which they are encountered
}
  \item{verbose}{
    If \code{TRUE}, this turns on output from the function. Otherwise
    only error and warning messages are produced
}
}
\details{
Each line on the input file represents a single call and is split into
fields using the function \code{strsplit}. The required fields are
extracted according to the \code{fields} argument. This \emph{must}
contain the locations of the sample and snp identifier
fields and \emph{either} the location of a genotype field \emph{or} the
locations of two allele fields.

If the \code{samples} and \code{snps} arguments contain vectors of
character strings, a \code{SnpMatrix} is created with these row and
column names and the  genotype values are "cherry-picked" from the input
file. If either, or both, of these arguments are specified simply as
numbers, then these
numbers determine the \emph{dimensions} of the \code{SnpMatrix}
created. In this case samples and/or SNPs are included in the
\code{SnpMatrix} on a first-come-first-served basis. If either
or both of these arguments are omitted, a preliminary scan of the input file
is carried out to find the missing sample and/or SNP identifiers. 
In this scan, 
when a sample or SNP identifier differs from that in the previous
line, but is identical to one previously found, then all the relevant
identifiers are assumed to have been found. This implies that
the file must be sorted, in some consistent order,
by sample and by SNP (although either one of these may vary fastest).

If the genotype is read as a single field, it's handling is specified
by the \code{gcodes} argument. If this is absent, \code{NULL}, or
\code{NA},  then
the genotype is assumed to be represented by a two-character field
(the two characters representing the two alleles). If \code{gcodes} is
a single string, then this is assumed to
be a regular expression which will split the genotype field into two allele
fields. Otherwise, \code{gcode} must be an array of length three,
specifying the three genotype codes in the order "AA", "AB", "BB"
}
\value{
If the genotype is read as a single field matching one of three
specified codes, the function returns an object of class
\code{SnpMatrix}. Otherwise it returns a list whose first element is the
\code{SnpMatrix} object and whose second element is a dataframe
containing the allele codes, with the SNP identifiers as row names. Note
that allele codes only occur in this file if they occur in a genotype
which was accepted. Thus, monomorphic SNPs have \code{allele.B} coded as
\code{NA}, and SNPs which never pass confidence score filters have both
alleles coded as \code{NA}.
}
\author{
David Clayton \email{david.clayton@cimr.cam.ac.uk}
}
\note{
  Unlike \code{read.snps.long},
  this function is written entirely in R and may not be particularly
  fast. However, it imposes no restrictions on the allele codes
  recognized.

  Homozygous genotypes are assumed to be represented in the input file 
  by coding both alleles to the same value. No special provision is made
  to read \code{XSnpMatrix} 
  objects; such data should first be read as a \code{SnpMatrix} and then
  coerced to an \code{XSnpMatrix} using \code{new} or \code{as}.
}
\seealso{
\code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}
}
\examples{
##
## No example supplied yet
##
}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
