\name{plotUncertainty}
\alias{plotUncertainty}
\title{
Plot posterior probabilities of genotype assignment
}
\description{
The snpStats package allows for storage of uncertain genotype
assignments in a one byte "raw" variable. The probabilities of
assignment form a three-vector, subject to the linear constraint that
they sum to 1.0; their possible values are grouped into 253 different
classes. This function displays counts of these classes on a
two-dimensional isometric plot.
}
\usage{
plotUncertainty(snp, nlevels = 10, color.palette = heat.colors(nlevels))
}
\arguments{
  \item{snp}{
    One or more columns of a \code{SnpMatrix} object
}
  \item{nlevels}{
    Probability cells are coloured according to frequency. This argument
    gives the number of colours that can be used
}
  \item{color.palette}{
    The colour palette to be used
}
}
\details{
The plot takes the form of an equilateral triangle in which each apex
represents a certain assignment to one of the three genotypes. A point
within the triangle represents, by the perpendicular distance from each
side, the three probabilities. Each of the 253 probability classes is
represented by a hexagonal cell, coloured according to its frequency in
the data, which is also written within the cell
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\examples{
## No example available yet
}
\keyword{hplot}

