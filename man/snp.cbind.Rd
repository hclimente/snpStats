\name{snp.cbind}
\alias{cbind}
\alias{rbind}
\alias{snp.cbind}
\alias{snp.rbind}
\alias{cbind2}
\alias{rbind2}
\title{Bind together two or more SnpMatrix objects}
\description{
These functions bind together two or more objects of class
\code{"SnpMatrix"} or \code{"XSnpMatrix"}. 
}
\usage{
# cbind(...)
# rbind(...)
snp.cbind(...)
snp.rbind(...)
}
\arguments{
  \item{\dots}{Objects of class \code{"SnpMatrix"} or \code{"XSnpMatrix"}.}
}
\details{
These functions  reproduce the action of the standard functions
\code{cbind} and \code{rbind}. These are constrained to
work by recursive calls to the generic functions  \code{cbind2}
and \code{rbind2} which take just two arguments. This
is somewhat inefficient in both time and memory use when binding more
than two objects, so the functions \code{snp.cbind} and
\code{snp.rbind}, which take multiple arguments, are also supplied.
  
When matrices are bound together by column, row names must be identical,
column names must not be duplicated and, for objects of class
\code{XSnpMatrix} the contents of the \code{Female} slot much match.
When matrices are bound by row,
column names must be identical and duplications of row names generate
warnings. 
}
\value{
  A new matrix, of the same type as the input matrices.
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link[base]{cbind}}, \code{\link[base:cbind]{rbind}}}
\examples{
data(testdata)
# subsetting ( Autosomes[c(1:9,11:19,21:29),] ) is quicker. this is just for illustrating
# rbind and cbind
first <- Autosomes[1:9,]
second <- Autosomes[11:19,]
third <- Autosomes[21:29,]
result1 <- rbind(first, second, third)
result2 <- snp.rbind(first, second, third)
all.equal(result1, result2)

result3 <- Autosomes[c(1:9,11:19,21:29),]
all.equal(result1, result3)

first <- Autosomes[,1:9]
second <- Autosomes[,11:19]
third <- Autosomes[,21:29]
result1 <- cbind(first, second, third)
result2 <- snp.cbind(first, second, third)
all.equal(result1, result2)

result3 <- Autosomes[,c(1:9,11:19,21:29)]
all.equal(result1, result3)

first <- Xchromosome[1:9,]
second <- Xchromosome[11:19,]
third <- Xchromosome[21:29,]
result1 <- rbind(first, second, third)
result2 <- snp.rbind(first, second, third)
all.equal(result1, result2)

result3 <- Xchromosome[c(1:9,11:19,21:29),]
all.equal(result1, result3)

first <- Xchromosome[,1:9]
second <- Xchromosome[,11:19]
third <- Xchromosome[,21:29]
result1 <- cbind(first, second, third)
result2 <- snp.cbind(first, second, third)
all.equal(result1, result2)

result3 <- Xchromosome[,c(1:9,11:19,21:29)]
all.equal(result1, result3)
}
\keyword{utilities}
