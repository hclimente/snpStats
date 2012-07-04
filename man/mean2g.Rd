\name{mean2g}
\alias{mean2g}
\alias{post2g}
\alias{g2post}
\title{
Raw coding of posterior probabilities of SNP genotype
}
\description{
An uncertain SNP genotype call is represented by the three posterior
probabilities of the three possible calls. In the class
\code{SnpMatrix}, an approximation to these is packed into a single
1-byte variable of type \code{raw}. These functions carry out this coding (and
decoding). 
}
\usage{
post2g(p, transpose = FALSE)
mean2g(m, maxE = FALSE)
g2post(g, transpose = FALSE)
}
\arguments{
  \item{p}{
    A matrix of posterior probabilities. If \code{transpose} is
    \code{FALSE} this is \code{Nx3}, otherwise \code{3xN} 
  }
  \item{m}{
    A vector of posterior means
  }
  \item{g}{
    A \code{raw} vector of genotype codes
  }
  \item{transpose}{
    A logical flag indicating transposition of the matric of posterior
    probabilities (see Description)
  }
  \item{maxE}{
    A logical flag selecting the maximum entropy option in
    \code{mean2g}
  }
}
\details{
\code{post2g} and \code{g2post} convert from posterior probabilities to
\code{raw} code and back respectively. If only the posterior expectation
of the genotype (when numerically coded 0, 1, or 2) is available, no
unique soultion exists in general and the behaviour of the function is
determined by the value of \code{maxE}. If \code{TRUE}, then the maximum
entropy solutions are returned while, if \code{FALSE}, an attempt is made
to return the least uncertain solution, by setting the posterior
probability of the \code{BB} genotype to zero when the posterior mean is
less than 1.0 and the probability of \code{AA} to zero when the mean is
greater than 1.0.
}
\value{
\code{post2g} and \code{mean2g} return a vector of type
\code{raw}. \code{g2post} returns a numeric matrix. 
}
\author{
David Clayton \email{dc208@cam.ac.uk}
}
\note{
These functions are provided mainly for users wishing to write their own
data input functions.
}
\examples{
data(testdata)
g <- Autosomes[1:10, 20] ## A vector of codes
p <- g2post(g)           ## Transform to probabilities ...
pg <- post2g(p)          ## ... and back to codes
m <- p[,2]+2*p[,3]       ## Posterior expectations
mg <- mean2g(m)          ## Transform to codes ...
pmg <- g2post(mg)        ## ... and transform to probabilities
## Write everything out
print(cbind(as.numeric(g), p, as.numeric(pg), m, as.numeric(mg), pmg))
}
\keyword{programming}

