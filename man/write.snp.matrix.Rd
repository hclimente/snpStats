\name{write.SnpMatrix}
\alias{write.SnpMatrix}
\title{Write a SnpMatrix object as a text file}
\description{
  This function is closely modelled on \code{write.table}. It writes an
  object of class \code{SnpMatrix} as a text file with one line for
  each row of the matrix. Genotpyes are written in numerical form,
  \emph{i.e.} as 0, 1 or 2 (where 1 denotes heterozygous) or,
  optionally, as a pair of alleles (each coded 1 or 2).
}
\usage{
write.SnpMatrix(x, file, as.alleles= FALSE, append = FALSE, quote = TRUE, sep = " ", eol = "\n", na = "NA", row.names = TRUE, col.names = TRUE)
}
\arguments{
  \item{x}{The object to be written }
  \item{file}{The name of the output file}
  \item{as.alleles}{If \code{TRUE}, write each genotype as two alleles}
  \item{append}{If \code{TRUE}, the output is appended to the designated
  file. Otherwise a new file is opened}
  \item{quote}{If \code{TRUE}, row and column names will be enclosed in quotes}
  \item{sep}{The string separating entries within a line}
  \item{eol}{The string terminating each line}
  \item{na}{The string written for missing genotypes}
  \item{row.names}{If \code{TRUE}, each row will commence with the row name}
  \item{col.names}{If \code{TRUE}, the first line will contain all the
    column names}
}
\value{
  A numeric vector giving the dimensions of the matrix written
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link[utils]{write.table}}, \code{\link{SnpMatrix-class}},
  \code{\link{XSnpMatrix-class}}}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
