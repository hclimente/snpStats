
\name{SnpMatrix-class}
\docType{class}
\alias{SnpMatrix-class}
\alias{initialize,SnpMatrix-method} 
\alias{[,SnpMatrix,ANY,ANY,missing-method}
\alias{rbind2,SnpMatrix,SnpMatrix-method}
\alias{cbind2,SnpMatrix,SnpMatrix-method}
\alias{coerce,SnpMatrix,numeric-method}
\alias{coerce,SnpMatrix,character-method}
\alias{coerce,matrix,SnpMatrix-method}
\alias{is.na,SnpMatrix-method}
\alias{cbind,SnpMatrix-method}
\alias{rbind,SnpMatrix-method}
\alias{show,SnpMatrix-method}
\alias{summary,SnpMatrix-method}
\alias{switch.alleles,SnpMatrix,ANY-method}


\title{Class "SnpMatrix"}
\description{ This class defines objects holding large arrays of single
  nucleotide polymorphism (SNP) genotypes generated using array
  technologies.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SnpMatrix", x)}
  where \code{x} is a matrix with storage mode \code{"raw"}.
	 Chips (usually corresponding to samples or
	 subjects) define rows of the matrix while polymorphisms (loci) define
	 columns. Rows and columns will usually have names which can be
	 used to link the data to further data concerning samples and SNPs
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} and storage
      mode \code{raw}  
    Internally, missing data are coded 0 and SNP genotypes are coded
    1, 2 or 3. Imputed values may not be known exactly. Such
    uncertain calls are grouped by probability and represented by codes
    4 to 253}
  }
}
\section{Extends}{
Class \code{"matrix"}, from data part.
Class \code{"structure"}, by class \code{"matrix"}.
Class \code{"array"}, by class \code{"matrix"}.
Class \code{"vector"}, by class "matrix", with explicit coerce.
Class \code{"vector"}, by class "matrix", with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x = "SnpMatrix", i = "ANY", j = "ANY",
	drop = "missing")}: subset
      operations}
    \item{cbind2}{\code{signature(x = "SnpMatrix", y = "SnpMatrix")}:
      S4 generic function to provide cbind() for two or
      more matrices together by column. Row names must match and column
      names must not coincide. If the matrices are of the derived class
      \code{\link{XSnpMatrix-class}}, the \code{Female} slot values must also
      agree}
    \item{coerce}{\code{signature(from = "SnpMatrix", to = "numeric")}:
      map to numeric values 0, 1, 2 or, for uncertain assignments, to
      the posterior expectation of the 0, 1, 2 code}
    \item{coerce}{\code{signature(from = "SnpMatrix", to =
	"character")}: map to codes "A/A", "A/B", "B/B", ""}
    \item{coerce}{\code{signature(from = "matrix", to = "SnpMatrix")}:
      maps numeric matrix (coded 0, 1, 2 or NA) to a SnpMatrix}
    \item{coerce}{\code{signature(from = "SnpMatrix", to =
	"XSnpMatrix")}:
      maps a SnpMatrix to an XSnpMatrix. Sex is inferred from the
      genotype data since males should not be heterozygous at any locus.
      After inferring sex, heterozygous calls for males are set to
      \code{NA}}
    \item{is.na}{\code{signature(x = "SnpMatrix")}: returns a logical
      matrix indicating whether each element is \code{NA}}
    \item{rbind2}{\code{signature(x = "SnpMatrix", y = "snp.matrix")}:
      S4 generic function to provide rbind() for two or
      more matrices by row. Column names must match and duplicated row
      names prompt warnings}
    \item{show}{\code{signature(object = "SnpMatrix")}: shows the size
      of the matrix (since most objects will be too large to show in full)}
    \item{summary}{\code{signature(object = "SnpMatrix")}: returns
      summaries of the data frames returned by
      \code{\link{row.summary}} and \code{\link{col.summary}}}
    \item{is.na}{\code{signature(x = "SnpMatrix")}: returns a logical
      matrix of missing call indicators}
    \item{switch.alleles}{\code{signature(x = "SnpMatrix", snps
	="ANY")}: Recode specified columns of of the matrix to reflect
      allele switches}

	 }
}
\references{\url{http://www-gene.cimr.cam.ac.uk/clayton}}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{This class requires at least version 2.3 of R}
\seealso{\code{\link{XSnpMatrix-class}}
}
\examples{
data(testdata)
summary(Autosomes)

# Just making it up - 3-10 will be made into NA during conversion
snps.class<-new("SnpMatrix", matrix(1:10))
snps.class
if(!isS4(snps.class)) stop("constructor is not working")

pretend.X <- as(Autosomes, 'XSnpMatrix')
if(!isS4(pretend.X)) stop("coersion to derived class is not S4")
if(class(pretend.X) != 'XSnpMatrix') stop("coersion to derived class is not working")

pretend.A <- as(Xchromosome, 'SnpMatrix')
if(!isS4(pretend.A)) stop("coersion to base class is not S4")
if(class(pretend.A) != 'SnpMatrix') stop("coersion to base class is not working")

# display the first 10 snps of the first 10 samples
print(as(Autosomes[1:10,1:10], 'character'))

# convert the empty strings (no-calls) explicitly to "NC" before
# writing to an (anonymous and temporary) csv file
csvfile <- tempfile()
write.csv(file=csvfile, gsub ('^$', 'NC',
                              as(Autosomes[1:10,1:10], 'character')
                             ), quote=FALSE)
unlink(csvfile)
}
\keyword{classes}
