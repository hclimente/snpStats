\name{XSnpMatrix-class}
\docType{class}
\alias{XSnpMatrix-class}
\alias{initialize,XSnpMatrix-method}
\alias{[,XSnpMatrix,ANY,ANY,ANY-method}
\alias{[<-,XSnpMatrix,ANY,ANY,XSnpMatrix-method}
\alias{coerce,XSnpMatrix,character-method}
\alias{coerce,SnpMatrix,XSnpMatrix-method}
\alias{show,XSnpMatrix-method}
\alias{summary,XSnpMatrix-method}

\title{Class "XSnpMatrix"}
\description{This class extends the  \code{\link{SnpMatrix-class}} to
    deal with SNPs on the X and Y chromosomes and mitocondrial SNPs.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("XSnpMatrix", x,
  diploid)}.
Such objects have an additional \code{slot} to objects of class
\code{"SnpMatrix"}
consisting of  a logical array of the same length as the number of
rows. This array indicates whether genotypes in that row are diploid
(\code{TRUE}) or haploid (\code{FALSE} as, for example,  SNPs on the X
chromosome for males).
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} and storage mode
      \code{"raw"}} 
    \item{\code{diploid}:}{Object of class \code{"logical"} indicating
      sex of samples}
  }
}
\section{Extends}{
Class \code{"SnpMatrix"}, directly, with explicit coerce.
Class \code{"matrix"}, by class \code{"SnpMatrix"}.
Class \code{"structure"}, by class \code{"SnpMatrix"}.
Class \code{"array"}, by class \code{"SnpMatrix"}.
Class \code{"vector"}, by class "SnpMatrix", with explicit coerce.
Class \code{"vector"}, by class "SnpMatrix", with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x = "XSnpMatrix", i = "ANY", j = "ANY",
	drop = "missing")}: subset extraction}
    \item{[<-}{\code{signature(x = "XSnpMatrix", i = "ANY", j = "ANY", "XSnpMatrix")}: subset
      assignment operation to replace part of an object}
    \item{coerce}{\code{signature(from = "XSnpMatrix", to =
	"character")}:  map to codes 0, 1, 2, or NA}
    \item{coerce}{\code{signature(from = "SnpMatrix", to =
	"XSnpMatrix")}:
      maps a SnpMatrix to an XSnpMatrix. Ploidy is inferred from the
      genotype data since haploid genotypes should always be coded as
      homozygous. After inferring ploidy, heterozygous calls for haploid
      genotpes are set to \code{NA}}
    \item{show}{\code{signature(object = "XSnpMatrix")}: map to codes
      "A/A", "A/B", "B/B", "A", "B" or ""}
    \item{summary}{\code{signature(object = "XSnpMatrix")}: returns
      the distribution of ploidy, together with
      summaries of the data frames returned by
      \code{\link{row.summary}} and \code{\link{col.summary}}}
  }
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{SnpMatrix-class}}
}
\examples{
data(testdata)
summary(Xchromosome)

# display the first 10 snps of the first 10 samples
print(as(Xchromosome[1:10,1:10],'character'))

# convert the empty strings (no-calls) explicitly to "NC" before
# writing to an (anonymous and temporary) csv file
csvfile <- tempfile()
write.csv(file=csvfile, gsub ('^$', 'NC',
                              as(Xchromosome[1:10,1:10], 'character')
                             ), quote=FALSE)
unlink(csvfile)
}
\keyword{classes}
