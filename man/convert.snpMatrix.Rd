\name{convert.snpMatrix}
\alias{convert.snpMatrix}
\alias{convert.snpMatrix.dir}
\title{
Convert \code{snpMatrix} objects to \code{snpAssoc} objects
}
\description{
These functions convert \code{snpMatrix} objects to \code{snpAssoc}
objects. \code{convert.snpMatrix} converts a single object, while
\code{convert.snpMatrix.dir} converts all stored elements in a specified
directory. They really only change the class names since most of the
classes in \code{snpAssoc} are backwards-compatible with 
\code{snpMatrix}. The exception is the \code{ImputationRules} class; 
\code{imputation.rules} objects will need to be regenerated.
}
\usage{
convert.snpMatrix(object)

convert.snpMatrix.dir(dir = ".", ext = "RData")
}
\arguments{
  \item{object}{
    Object to be converted
  }
  \item{dir}{
    A directory containing saved \code{snpMatrix} objects
  }
  \item{ext}{
    The file extension for files containing such objects
  }
}
\value{
\code{convert.snpMatrix} returns the converted object. 
\code{convert.snpMatrix.dir} rewrites the files in place.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\keyword{classes}
