\name{write.plink}
\alias{write.plink}
\title{
Write files for analysis in the PLINK toolset
}
\description{
Given a \code{SnpMatrix} object, together with associated subject and
SNP support dataframes, this function writes \code{.bed}, 
\code{.bim}, and \code{.fam} files for processing in the PLINK toolset
}
\usage{
write.plink(file.base, snp.major = TRUE, snps,
   subject.data, pedigree, id, father, mother, sex, phenotype,
   snp.data, chromosome, genetic.distance, position, allele.1, allele.2,
   na.code = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.base}{A character string giving the base filename. The
    extensions \code{.bed}, \code{.bim}, and \code{.fam} are appended to
    this string to give the filenames of the three output files
  }
  \item{snp.major}{Logical variable controlling whether the \code{.bed}
    file is in SNP-major or subject-major order 
  }
  \item{snps}{The \code{SnpMatrix} or \code{XSnpMatrix} object to be written out
  }
  \item{subject.data}{(Optional) A subject support dataframe. If
    supplied, the next six arguments (which define the fields of the
    PLINK \code{.fam} file) will be evaluated in this
    environment, after matching row names with the row names of \code{snps}.
    Otherwise they will be evaluated in the calling environment; they then
    must be of the right length and in the correct order.
  }
  \item{pedigree}{A pedigree (family) identifier. Default is the row
    names of \code{snps}.
  }
  \item{id}{An identifier of an individual within family. Default is a
    vector of \code{na.code}.
  }
  \item{father}{The within-family identifier of the subject's
    father. Default is a vector of  \code{na.code}.
  }
  \item{mother}{The within-family identifier of the subject's
    mother. Default is a vector of  \code{na.code}.
  }
  \item{sex}{Sex of the individual. Default is a vector of
    \code{na.code}. This will be coerced to type \code{numeric}.
  }
  \item{phenotype}{The primary phenotype value. Default is a vector of
     \code{na.code}. This will be coerced to type \code{numeric}.
  }
  \item{snp.data}{(Optional) A SNP support dataframe. If
    supplied, the next five arguments (which define the  columns of the
    PLINK \code{.bim} file) will be evaluated in this
    environment,  after matching row names with the column names of
    \code{snps}. Otherwise they will be evaluated in the calling
    environment;  they then
    must be of the right length and in the correct order. 
  }
  \item{chromosome}{The chromosome on which the SNP is located. This
    should either be numeric, as specified by SPLINK, or character, with
    \code{"X"}, \code{"Y"}, \code{"XY"}, and \code{"MT"} for the
    non-numeric values. Default 
    is a vector of  \code{na.code}, or a vector of 23's if \code{snps}
    is a \code{XSnpMatrix}.
  }
  \item{genetic.distance}{The location of the SNP, expressed as a
    genetic distance. Default is a vector of  \code{na.code}.
    This will be coerced to type \code{numeric}.
  }
  \item{position}{The physical location of the SNP, expressed in base pairs.
    Default is a vector of  \code{na.code}.
    This will be coerced to type \code{numeric}.
  }
  \item{allele.1}{A character vector giving the first allele. Default is
    a vector of \code{"A"}s.
  }
  \item{allele.2}{A character vector giving the first allele. Default is
    a vector of \code{"B"}s.
  }
  \item{na.code}{The code to be written for \code{NA}
    in the \code{.fam} and \code{.bin} output files. It should be
    numeric (or capable of coercion to numeric).
  }
}
\details{
  For more details of required codings in \code{.fam} and \code{.bim}
  files, see the PLINK documentation.
}
\value{
Returns \code{NULL}.
}
\references{
  PLINK: Whole genome association analysis toolset.
  \url{http://pngu.mgh.harvard.edu/~purcell/plink/}
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{read.plink}},
  \code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}
}
\examples{
data(testdata)
## the use of as.numeric() below is not necessary since this is done
## automatically. It just illustrates use of expressions for these arguments
## Note that cc and sex are variables within the subject.data frame
## This command writes files test.bed. test.fam and test.bim
write.plink(file.base="test", snps=Autosomes,
    subject.data=subject.data, phenotype = as.numeric(cc), sex=as.numeric(sex),
    snp.major=FALSE)
}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}

