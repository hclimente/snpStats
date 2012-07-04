\name{tdt.snp}
\alias{tdt.snp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{1-df and 2-df tests for genetic associations with SNPs (or
  imputed SNPs) in family data}
\description{
Given large-scale SNP data for families comprising both parents and one or 
more affected offspring, this function computes 1 df tests (the TDT test) and 
a 2 df test based on observed and expected transmissions of genotypes. 
Tests based on imputation rules can also be carried out.
}
\usage{
tdt.snp(ped, id, father, mother, affected, data = sys.parent(), snp.data, 
    rules = NULL, snp.subset, check.inheritance = TRUE, robust = FALSE,
    uncertain = FALSE, score = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ped}{Pedigree identifiers}
  \item{id}{Subject identifiers}
  \item{father}{Identifiers for subjects' fathers}
  \item{mother}{Identifiers for subjects' mothers}
  \item{affected}{Disease status (TRUE if affected, FALSE otherwise)}
  \item{data}{A data frame in which to evaluate the previous five arguments}
  \item{snp.data}{An object of class \code{"SnpMatrix"} containing the SNP
    genotypes to be tested}
  \item{rules}{An object of class
    \code{"ImputationRules"}. If
    supplied, the rules coded in this object are used, together with
    \code{snp.data}, to calculate tests for imputed SNPs}
  \item{snp.subset}{A vector describing the subset of SNPs to be
    considered. Default action is to test all SNPs in \code{snp.data}
    or, in imputation mode,  as specified by \code{rules}}
  \item{check.inheritance}{If TRUE, each affected offspring/parent trio 
    is tested for Mendelian inheritance and excluded if the test fails. If
    FALSE, misinheriting trios are used but the "robust" variance option is
    forced}
  \item{robust}{If TRUE, forces the robust (Huber-White) variance option 
    (with \code{ped} determining independent "clusters")}
  \item{uncertain}{If TRUE, uncertain genotypes are handed by replacing
    score contributions by their posterior expectations. Otherwise these
    are treated as missing. Setting this option authomatically invokes
    use of \code{robust} variance estimates}
  \item{score}{If \code{TRUE}, the output object will contain, for each
    SNP,  the score vector and its variance-covariance matrix}
}
\details{
  Formally, the test statistics are score tests for the "conditioning on
  parental genotype" (CPG) likelihood. Parametrization of associations is
  the same as for the population-based tests calculated by 
  \code{\link{single.snp.tests}} so that results from family-based and 
  population-based studies can be combined using \code{\link{pool}}.

  When the function is used to calculate tests for imputed SNPs, the
  test is still an approximate score test. The current version does not use 
  the family relationships in the imputation. With this option, the robust 
  variance estimate is forced. 
  
  The first five arguments are usually derived from  a "pedfile". If a 
  data frame is supplied for the \code{data} argument, the first five 
  arguments will be evaluated in this frame. Otherwise they will be evaluated 
  in the calling environment. If the arguments are missing, they will be 
  assumed to be in their usual positions in the pedfile data frame 
  i.e. in columns one to 
  four for the identifiers  and column six for disease status 
  (with affected coded \code{2}). If the pedfile data are obtained from
  a dataframe, the row names of the \code{data} and \code{snp.data}
  files will be used to align the pedfile and SNP data. Otherwise, these
  vectors will be assumed to be in the same order as the rows of
  \code{snp.data}. 

  The \code{snp.subset} argument can be a logical,
  integer, or character vector.

  If imputed rather than observed SNPs are tested, or
  if \code{check.inheritance} is set to
  \code{FALSE}, the "robust" variance estimate is used regardless of the
  value supplied for the \code{robust} argument.
}
\value{
  An object of class
  \code{"SingleSnpTests"}. 
  If \code{score=TRUE}, the output object will be of the extended class
  \code{"SingleSnpTestsScore"}
  containing additional slots holding the score statistics and their
  variances (and covariances). This allows meta-analysis using the
  \code{\link{pool}} function.
}
\references{
  Clayton (2008) Testing for association on the X chromosome
  \emph{Biostatistics}, \bold{9}:593-600.)
}
\note{When the snps are on the X chromosome (i.e. when the \code{snp.data} 
  argument is of class \code{"XSnpMatrix"}), the tests are constructed 
  in the same way as was described by Clayton (2008) for population-based 
  association tests i.e. assuming that 
  genotype relative risks for males mirror thos of homozygous females
}
\author{David Clayton \email{dc208@cam.ac.uk}}
\seealso{\code{\link{single.snp.tests}}, \code{\link{impute.snps}}, 
  \code{\link{pool}}, \code{\link{ImputationRules-class}},
  \code{\link{SingleSnpTests-class}},
  \code{\link{SingleSnpTestsScore-class}} 
}
\examples{
data(families)
tdt.snp(data=pedData, snp.data=genotypes)
}
\keyword{htest}
