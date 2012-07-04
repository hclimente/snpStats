\name{test.allele.switch}
\alias{test.allele.switch}
\title{Test for switch of alleles between two collections}
\description{
When testing genotype data derived from different platforms or scoring
algorithms a common problem is switching of alleles. This function
provides a diagnostic for this. Input can either be two objects of class
\code{"SnpMatrix"} to be examined, column by column, for allele
switching, or a single \code{"SnpMatrix"} object together with an
indicator vector giving group membership for its rows.
}
\usage{
test.allele.switch(snps, snps2 = NULL, split = NULL, prior.df = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{An object of class \code{"SnpMatrix"} or \code{"XSnpMatrix"}}
  \item{snps2}{A second object of the same class as \code{snps}}
  \item{split}{If only one \code{SnpMatrix} object supplied, a vector
    with the same number of elements as rows of \code{snps}. It must be
    capable of coercion to a factor with two levels.
  }
  \item{prior.df}{A degree of freedom parameter for the prior
    distribution of the allele frequency \code{prior.df} (see Details)}
}
\details{
This function calculates a Bayes factor for the comparison of the
hypothesis that the alleles have been switched with the hypothesis that
they have not been switched. This requires integration over the
posterior distribution of the allele frequency. The prior is taken as a
beta distribution with both parameters equal to \code{prior.df}so that
the prior is symmetric about 0.5. The default,
\code{prior.df=1} represents a uniform prior on (0,1). 
}
\value{
A vector containing the log (base 10) of the Bayes Factors for an allele
switch. 
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{SnpMatrix-class}}, \code{\link{XSnpMatrix-class}}}
\examples{
data(testdata)
#
# Call with two SnpMatrix arguments
#
cc <- as.numeric(subject.data$cc)
lbf1 <- test.allele.switch(Autosomes[cc==1,], Autosomes[cc==2,])
#
# Single matrix call (giving the same result)
#
lbf2 <- test.allele.switch(Autosomes, split=cc)
}
\keyword{utilities}
